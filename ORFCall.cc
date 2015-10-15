///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ORFCall.cc
 *
 *  Created on: Feb 27, 2014
 *      Author: tsharpe
 */
#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstring>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "gzstream.h"

#define FatalErr(message) std::cerr<<message<<std::endl,exit(1)

namespace
{
char const FASTQ_QUAL_OFFSET = 33;
unsigned const NCODONS = 64;

typedef unsigned char Call;
inline Call complement( Call call )
{ return (call==4) ? 4 : (call^3); }
inline char charFor( Call call )
{ return "ACGTN"[call]; }

typedef std::vector<Call> Sequence;
typedef Sequence::size_type Offset;
typedef unsigned char Qual;
typedef std::vector<Qual> Quals;

// describes an exon as a pair of offsets on the reference sequence.
struct Exon
{
    explicit Exon( Offset beg, Offset end=0 )
    : mBegin(beg), mEnd(end) {}

    Offset size() const { return mEnd-mBegin; }

    Offset mBegin;
    Offset mEnd;
};
typedef std::vector<Exon> ORF;

typedef Sequence::const_iterator Kmer;
// produces a hash value for the next K bases from a sequence.
class KmerHasher
{
public:
    KmerHasher( unsigned K ) : mK(K) {}

    size_t operator()( Kmer itr ) const
    { size_t result = 14695981039346656037ul;
      for ( Kmer end(itr+mK); itr != end; ++itr )
          result = 1099511628211ul*(result ^ *itr); // FNV-1a algorithm
      return result; }

private:
    unsigned const mK;
};
// compares K bases from 2 sequences for equality.
class KmerComparator
{
public:
    KmerComparator( unsigned K ) : mK(K) {}

    bool operator()( Kmer itr1, Kmer itr2 ) const
    { return std::equal(itr1,itr1+mK,itr2); }

private:
    unsigned const mK;
};
// describes a kmer as unique within the reference or duplicated,
// and, if unique, where it occurs within the reference and its orientation.
class KmerLoc
{
    enum class Status : char
    { UNSET, UNIQUE, DUP };

public:
    KmerLoc() : mStatus(Status::UNSET) {}

    Offset getOffset() const { return mOffset; }
    bool isRC() const { return mIsRC; }
    bool isUnique() const { return mStatus == Status::UNIQUE; }

    void set( unsigned offset, bool isRC )
    { if ( mStatus != Status::UNSET ) mStatus = Status::DUP;
      else { mOffset = offset; mIsRC = isRC; mStatus = Status::UNIQUE; } }

private:
    Offset mOffset;
    bool mIsRC;
    Status mStatus;
};
// a map of kmers onto KmerLocs
typedef std::unordered_map<Kmer,KmerLoc,KmerHasher,KmerComparator> KmerDict;

// describes a gap-free alignment of a sequence onto the reference
class Alignment
{
public:
    Alignment() : mRefStart(0), mRefEnd(0), mReadStart(0), mReadRC(false) {}
    Alignment( Offset refStart, Offset refEnd, Offset readStart, bool readRC )
    : mRefStart(refStart), mRefEnd(refEnd),
      mReadStart(readStart), mReadRC(readRC) {}

    Alignment& operator+=( Offset off )
    { mRefStart += off; mReadStart += off; return *this; }

    bool empty() const { return mRefStart==mRefEnd; }
    Offset size() const { return mRefEnd-mRefStart; }
    Offset getRefStart() const { return mRefStart; }
    Offset getRefEnd() const { return mRefEnd; }
    Offset getReadStart() const { return mReadStart; }
    Offset getReadEnd() const { return mReadStart+(mRefEnd-mRefStart); }
    bool isReadRC() const { return mReadRC; }

    static bool disjoint( Alignment const& aln1, Alignment const& aln2 )
    { return aln1.mRefStart > aln2.mRefEnd || aln2.mRefStart > aln1.mRefEnd; }

private:
    Offset mRefStart;
    Offset mRefEnd;
    Offset mReadStart;
    bool mReadRC;
};

// everything we need to know about the sequenced construct.
class Reference
{
public:
    Reference( char const* refFile,
                unsigned const K, Qual minQ, double maxMismatchFreq )
    : mK(K), mMinQ(minQ), mMaxMismatchFrac(maxMismatchFreq),
      mUniqueKmers(0,KmerHasher(K),KmerComparator(K))
    { readRefSeq(refFile);
      if ( mRefSeq.size() < K )
          FatalErr("There are only " << mRefSeq.size()
                      << "bases in the reference file");
      buildRC();
      buildDictionary(K); }

    unsigned getK() const { return mK; }
    Sequence const& getRef() const { return mRefSeq; }

    ORF const& getORF() const { return mORF; }
    size_t getORFSize() const
    { size_t result = 0;
      for ( Exon const& exon : mORF ) result += exon.size();
      return result; }
    Sequence getORFSeq() const
    { Sequence result; result.reserve(getORFSize());
      for ( Exon const& exon : mORF )
        std::copy(mRefSeq.begin()+exon.mBegin,mRefSeq.begin()+exon.mEnd,
                    std::back_inserter(result) );
      return result; }

    bool findFirstAlignment( Sequence const& seq, Alignment* pAlignment ) const;
    bool findLastAlignment( Sequence const& seq, Quals const& quals,
                            Alignment* pAlignment ) const;

    void report() const
    { size_t refLen = mRefSeq.size();
      std::cerr << "Reference length: " << refLen << ".\n";
      std::cerr << "Dict size: " << mUniqueKmers.size()
                        << " (" << 2*(refLen-mK+1)-mUniqueKmers.size()
                        << " duplicate kmers).\n";
      std::cerr << "ORF length: " << getORFSize() << ".\n";
      std::cerr << "ORF segments:";
      char prefix = ' ';
      // leave final sentinel out of description (that's the end()-1 business.)
      for ( auto itr=mORF.begin(),end=mORF.end()-1; itr != end; ++itr )
      { std::cerr << prefix << itr->mBegin+1 << ':' << itr->mEnd+1;
        prefix = ','; }
      std::cerr << ".\n"; }

private:
    void readRefSeq( char const* refFile );

    void buildRC()
    { // store the RC of the reference sequence
      mRefSeqRC.reserve(mRefSeq.size());
      for ( auto itr=mRefSeq.rbegin(),end=mRefSeq.rend(); itr != end; ++itr )
        mRefSeqRC.push_back(complement(*itr)); }

    void buildDictionary( unsigned const K );
    double mismatchFraction( Sequence::const_iterator itr,
                             Sequence::const_iterator end,
                             Sequence::const_iterator rItr,
                             Quals::const_iterator qItr ) const
    { unsigned comparisons = 0, mismatches = 0;
      while ( itr != end )
      { if ( *qItr >= mMinQ ) { ++comparisons; mismatches += (*itr != *rItr); }
        ++itr; ++qItr; ++rItr; }
      return comparisons ? 1.*mismatches/comparisons : 1.; }

    unsigned mK;
    Qual mMinQ;
    double const mMaxMismatchFrac;
    Sequence mRefSeq;
    Sequence mRefSeqRC;
    ORF mORF;
    KmerDict mUniqueKmers;
};

void Reference::readRefSeq( char const* refFile )
{
    std::ifstream ifs(refFile,std::ios_base::in|
                                std::ios_base::binary|
                                std::ios_base::ate);
    if ( !ifs )
        FatalErr("Can't open reference file: " << refFile);
    mRefSeq.reserve(ifs.tellg());
    ifs.seekg(0);
    std::string line;
    if ( !std::getline(ifs,line) )
        FatalErr("Can't read reference file.");
    if ( line.empty() || line[0] != '>' )
        FatalErr("Reference file doesn't begin with a fasta comment.");
    unsigned lineNo = 0;
    bool inORF = false;
    while ( std::getline(ifs,line) )
    {
        ++lineNo;
        for ( char chr : line )
        {
            if ( std::isspace(chr) )
                continue;

            switch (chr)
            {
            case 'a': case 'A': mRefSeq.push_back(0); break;
            case 'c': case 'C': mRefSeq.push_back(1); break;
            case 'g': case 'G': mRefSeq.push_back(2); break;
            case 't': case 'T':
            case 'u': case 'U': mRefSeq.push_back(3); break;
            case '[':
                if ( inORF )
                    FatalErr("Nested ORF segment at line " << lineNo);
                inORF=true; mORF.emplace_back(mRefSeq.size()); break;
            case ']':
                if ( !inORF )
                    FatalErr("Unopened ORF segment at line " << lineNo);
                inORF=false; mORF.back().mEnd = mRefSeq.size(); break;
            default:
                FatalErr("Gibberish in reference file line " << lineNo << ":\n"
                        << line);
            }
        }
    }
    if ( inORF )
        FatalErr("Unclosed ORF segment.");
    if ( mORF.empty() )
    {
        std::cerr << "No ORF marked -- using entire reference sequence as ORF"
                << std::endl;
        mORF.emplace_back(0,mRefSeq.size());
    }
    if ( !ifs.eof() )
        FatalErr("Failed to read reference file to the end.");
    // add sentinel
    mORF.emplace_back(mRefSeq.size(),mRefSeq.size());
}

bool Reference::findFirstAlignment( Sequence const& seq, Alignment* pAln ) const
{
    if ( seq.size() < mK )
        return false;

    Kmer kmer(seq.begin());
    Kmer end(seq.end()-mK+1);
    auto notFound = mUniqueKmers.end();
    for ( auto kmer=seq.begin(),end=seq.end()-mK+1; kmer != end; ++kmer )
    {
        auto itr = mUniqueKmers.find(kmer);
        if ( itr != notFound )
        {
            Offset readStart = kmer-seq.begin();
            Offset refStart = itr->second.getOffset();
            bool isRC = itr->second.isRC();
            if ( isRC )
            {
                readStart = seq.size()-(readStart+mK);
                refStart = mRefSeq.size()-(refStart+mK);
            }
            *pAln = Alignment(refStart,refStart+mK,readStart,isRC);
            return true;
        }
    }
    return false;
}

bool Reference::findLastAlignment( Sequence const& seq, Quals const& quals,
                                    Alignment* pAln ) const
{
    if ( seq.size() < mK )
    {
        *pAln = Alignment();
        return false;
    }
    Kmer kmer(seq.end()-mK+1);
    Kmer beg(seq.begin());
    auto notFound = mUniqueKmers.end();
    while ( kmer != beg )
    {
        auto itr = mUniqueKmers.find(--kmer);
        if ( itr != notFound )
        {
            Offset refStart = itr->second.getOffset();
            Offset refEnd = refStart + mK;
            Offset readStart = kmer-beg;
            Offset backOff = std::min(refStart,readStart);
            refStart -= backOff;
            readStart -= backOff;
            bool isRC = itr->second.isRC();
            Sequence const& ref = isRC ? mRefSeqRC : mRefSeq;
            double mmFrac = mismatchFraction(ref.begin()+refStart,
                                              ref.begin()+refEnd,
                                              seq.begin()+readStart,
                                              quals.begin()+readStart);
            if ( isRC )
            {
                readStart = seq.size()-(readStart+(refEnd-refStart));
                refStart = mRefSeq.size()-refStart;
                refEnd = mRefSeq.size()-refEnd;
                std::swap(refStart,refEnd);
            }
            *pAln = Alignment(refStart,refEnd,readStart,isRC);
            return mmFrac <= mMaxMismatchFrac;
        }
    }
    *pAln = Alignment();
    return false;
}

void Reference::buildDictionary( unsigned const K )
{
    // build a dictionary of unique kmers in ref and ref RC
    mUniqueKmers.reserve(2*mRefSeq.size());
    Kmer end(mRefSeq.end()-K+1);
    Offset offset = 0;
    for ( Kmer itr(mRefSeq.begin()); itr != end; ++itr )
        mUniqueKmers[itr].set(offset++,false);

    end = mRefSeqRC.end()-K+1;
    offset = 0;
    for ( Kmer itr(mRefSeqRC.begin()); itr != end; ++itr )
        mUniqueKmers[itr].set(offset++,true);

    // delete duplicates
    auto dItr(mUniqueKmers.begin()), dEnd(mUniqueKmers.end());
    while ( dItr != dEnd )
        if ( !dItr->second.isUnique() )
            dItr = mUniqueKmers.erase(dItr);
        else
            ++dItr;
}

// a parser for fastq-formatted files.
class FastqReader
{
public:
    FastqReader( char const* fastqFile )
    : mpIS(openStream(fastqFile)), mIS(*mpIS), mLineNo(0)
    { if ( !mIS ) FatalErr("Opening fastq file: " << fastqFile); }

    ~FastqReader()
    { delete mpIS; }

    bool getChunk( std::vector<Sequence>& seqs, std::vector<Quals>& quals,
                    unsigned chunkSize = 100000 );

private:
    std::istream* openStream( char const* fastqFile )
    { auto len = strlen(fastqFile);
      if ( len >= 3 && !memcmp(fastqFile+len-3,".gz",3) )
        return new igzstream(fastqFile,std::ios_base::in|std::ios_base::binary);
      return new std::ifstream(fastqFile,std::ios_base::in|std::ios_base::binary);
    }

    std::istream* mpIS;
    std::istream& mIS;
    unsigned mLineNo;
};

// fastq files are read in chunks to limit memory use
bool FastqReader::getChunk( std::vector<Sequence>& seqs,
                            std::vector<Quals>& quals,
                            unsigned chunkSize )
{
    if ( mIS.eof() )
        return false;

    seqs.clear(); seqs.reserve(chunkSize);
    quals.clear(); quals.reserve(chunkSize);

    std::string seqheader;
    std::string sSeq;
    std::string qualheader;
    std::string sQual;
    Sequence sq;
    Quals qs;
    while ( ++mLineNo && std::getline(mIS,seqheader) )
    {
        if ( seqheader.empty() )
            FatalErr("Read empty header at line " << mLineNo);
        if ( ++mLineNo && !std::getline(mIS,sSeq) )
            FatalErr("Early EOF at line " << mLineNo);
        if ( !sSeq.empty() && sSeq.back() == '\r' )
            sSeq.resize(sSeq.size()-1);
        if ( sSeq.empty() )
            FatalErr("Empty sequence at line " << mLineNo);
        if ( ++mLineNo && !std::getline(mIS,qualheader) )
            FatalErr("Early EOF at line " << mLineNo);
        if ( qualheader.empty() || qualheader[0] != '+' )
            FatalErr("Quals header doesn't begin with '+' at line " << mLineNo);
        if ( ++mLineNo && !std::getline(mIS,sQual) )
            FatalErr("Early EOF at line " << mLineNo);
        if ( !sQual.empty() && sQual.back() == '\r' )
            sQual.resize(sQual.size()-1);
        if ( sSeq.size() != sQual.size() )
            FatalErr("Calls and quals of unequal length at line " << mLineNo);

        sq.clear(); sq.reserve(sSeq.size());
        qs.clear(); qs.reserve(sSeq.size());
        auto qItr(sQual.begin());
        for ( char chr : sSeq )
        {
            switch (chr)
            {
            case 'a': case 'A': sq.push_back(0); break;
            case 'c': case 'C': sq.push_back(1); break;
            case 'g': case 'G': sq.push_back(2); break;
            case 't': case 'T':
            case 'u': case 'U': sq.push_back(3); break;
            default: sq.push_back(4); break;
            }
            qs.push_back(*qItr-FASTQ_QUAL_OFFSET);
            ++qItr;
        }
        seqs.push_back(sq);
        quals.push_back(qs);
        if ( seqs.size() == chunkSize )
            return true;
    }
    if ( !mIS.eof() )
        FatalErr("Failed to read to EOF (read to line " << mLineNo << ").");
    return !seqs.empty();
}

// a translator of codons (6-bit values encoding 3 base calls) into amino acids.
class AminoAcidNamer
{
public:
    AminoAcidNamer( std::string const& codonTranslation )
    : mLabels(codonTranslation)
    { if ( codonTranslation.size() != NCODONS )
        FatalErr("Codon translation string must have exactly " << NCODONS
                    << " characters.");
      auto beg = mLabels.begin(), end = mLabels.end();
      std::sort(beg,end);
      mLabels.resize(std::unique(beg,end)-beg);
      end = mLabels.end();
      for ( size_t idx = 0; idx != NCODONS; ++idx )
        mCodonToAA[idx] = std::lower_bound(beg,end,codonTranslation[idx])-beg; }

    void getAACounts( size_t const* codonCounts,
                            std::vector<size_t>* pAACounts ) const
    { std::vector<size_t>& aaCounts = *pAACounts;
      aaCounts.clear(); aaCounts.resize(mLabels.size());
      for ( size_t idx : mCodonToAA )
          aaCounts[idx] += *codonCounts++; }

    size_t indexFor( Sequence::const_iterator itr ) const
    { unsigned codon = *itr << 4; codon |= *++itr << 2; codon |= *++itr;
      return mCodonToAA[codon]; }

    // translate a 6-bit codon value to an amino acid index
    size_t indexFor( unsigned codon ) const { return mCodonToAA[codon]; }

    // a container of amino acid names
    size_t size() const { return mLabels.size(); }
    std::string::const_iterator begin() const { return mLabels.begin(); }
    std::string::const_iterator end() const { return mLabels.end(); }
    char operator[]( size_t idx ) const { return mLabels[idx]; }

private:
    std::string mLabels;
    size_t mCodonToAA[NCODONS];
};

struct CodonCounts
{
    CodonCounts() : mCounts{}, mIndelCounts{} {}

    size_t mCounts[NCODONS];
    size_t mIndelCounts[2];
};
typedef std::vector<CodonCounts> CodonCountsVec;
typedef CodonCountsVec::iterator CCVItr;

// keeps track of codon boundaries during sequential processing of sequence.
// bumps the appropriate counter when we have a whole codon.
class CodonAccumulator
{
public:
    CodonAccumulator( CCVItr beg, Offset orfOffset, Qual minQ )
    : mCCVItr(beg+orfOffset/3), mCount(orfOffset%3), mCodon(0), mIsHQ(!mCount),
      mMinQ(minQ) {}

    void operator()( Call call, Qual qual )
    { if ( (mIsHQ = mIsHQ && qual >= mMinQ) )
        mCodon = ((mCodon << 2) | call) & 0x3F;
      if ( ++mCount == 3 )
      { if ( mIsHQ ) mCCVItr->mCounts[mCodon] += 1;
        mIsHQ = true;
        mCount = 0;
        ++mCCVItr; } }

private:
    CCVItr mCCVItr;  // the counter to bump when the codon is complete
    unsigned mCount; // which position of the codon are we pointing to
    unsigned mCodon; // the 6-bit value of the codon
    bool mIsHQ;      // are all the calls in the codon of high quality?
    Qual mMinQ;      // what it means to be high quality
};

class SmithWatermanizer
{
public:
    // direction of traceback
    enum class Trace : unsigned char { DONE, UP, LEFT, DIAG };

    // the length of a base-per-base matchup, or of an indel
    struct Block
    {
        Block( Trace dir, unsigned len ) : mDir(dir), mLen(len) {}
        Trace mDir;
        unsigned mLen;
    };

    int align( Sequence::const_iterator bX, Sequence::const_iterator eX,
                Sequence::const_iterator bY, Sequence::const_iterator eY,
                std::ostream* pOS=nullptr )
    { int result = score(bX,eX,bY,eY);
      traceback(bX,eX,bY,eY,pOS);
      return result; }

    std::vector<Block> const& getMatchup() const
    { return mMatchup; }

private:
    static int const GAP_OPEN = -14;
    static int const GAP_EXTEND = -3;
    static int const MATCH = 4;
    static int const MIS_MATCH = -6;

    int score( Sequence::const_iterator bX, Sequence::const_iterator eX,
               Sequence::const_iterator bY, Sequence::const_iterator eY );

    void traceback( Sequence::const_iterator bX, Sequence::const_iterator eX,
                    Sequence::const_iterator bY, Sequence::const_iterator eY,
                    std::ostream* pOS );

    static int max3( int dScore, int hScore, int vScore, Trace* pTrace )
    { int result;
      if ( dScore >= hScore )
      { if ( dScore >= vScore )
        { *pTrace = Trace::DIAG; result = dScore; }
        else if ( hScore >= vScore )
        { *pTrace = Trace::LEFT; result = hScore; }
        else
        { *pTrace = Trace::UP; result = vScore; } }
      else if ( hScore >= vScore )
      { *pTrace = Trace::LEFT; result = hScore; }
      else
      { *pTrace = Trace::UP; result = vScore; }
      return result; }

    // row of values for a vertical move (introduces or extends a gap in X)
    std::vector<int> mMaxScores;

    // row of best scores considering all 3 moves for each column
    std::vector<int> mVScores;

    // traceback pointers for the entire matrix
    std::vector<Trace> mTraceback;

    // traceback builds this structure
    std::vector<Block> mMatchup;

    // buffers for printing alignments
    std::string mLineX;
    std::string mLineY;
    std::string mLineZ;
};
int const SmithWatermanizer::GAP_OPEN;
int const SmithWatermanizer::GAP_EXTEND;
int const SmithWatermanizer::MATCH;
int const SmithWatermanizer::MIS_MATCH;

int SmithWatermanizer::score( Sequence::const_iterator bX,
                              Sequence::const_iterator eX,
                              Sequence::const_iterator bY,
                              Sequence::const_iterator eY )
{
    using std::distance;
    auto nCols = distance(bX,eX);
    mVScores.resize(nCols);
    mMaxScores.resize(nCols);
    mTraceback.resize((nCols+1)*(distance(bY,eY)+1));
    Trace* pTrace = &mTraceback.front();
    *pTrace++ = Trace::DONE;
    std::fill(pTrace,pTrace+nCols,Trace::LEFT);
    pTrace += nCols;
    int rowInit = GAP_OPEN;
    int colInit = GAP_OPEN;
    Call baseX0 = *bX;

    Call baseY = *bY;
    int* pVScore = &mVScores.front();
    int* pMaxScore = &mMaxScores.front();
    *pTrace++ = Trace::UP;

    // (0,0) can't extend any gaps
    int dScore = (baseX0==baseY ? MATCH : MIS_MATCH);
    int hScore = rowInit + GAP_OPEN;
    int vScore = colInit + GAP_OPEN;
    int maxScore = max3(dScore,hScore,vScore,pTrace++);

    *pVScore++ = vScore;
    *pMaxScore++ = maxScore;

    // (1..n,0) can't extend vGap
    for ( Sequence::const_iterator iX=++bX; iX != eX; ++iX )
    {
        dScore = colInit + (*iX==baseY ? MATCH : MIS_MATCH);
        if ( pTrace[-1] == Trace::LEFT )
            hScore = std::max(maxScore+GAP_OPEN,hScore+GAP_EXTEND);
        else
            hScore = maxScore + GAP_OPEN;
        colInit += GAP_EXTEND;
        vScore = colInit + GAP_OPEN;
        maxScore = max3(dScore,hScore,vScore,pTrace++);

        *pVScore++ = vScore;
        *pMaxScore++ = maxScore;
    }

    int maxScorePrevRow;
    while ( ++bY != eY )
    {
        baseY = *bY;
        pVScore = &mVScores.front();
        pMaxScore = &mMaxScores.front();
        *pTrace++ = Trace::UP;

        // (0,1..m) can't extend hGap
        dScore = rowInit + (baseX0==baseY ? MATCH : MIS_MATCH);
        rowInit += GAP_EXTEND;
        hScore = rowInit + GAP_OPEN;
        maxScorePrevRow = *pMaxScore;
        vScore = std::max(maxScorePrevRow+GAP_OPEN,*pVScore+GAP_EXTEND);
        maxScore = max3(dScore,hScore,vScore,pTrace++);

        *pVScore++ = vScore;
        *pMaxScore++ = maxScore;

        // (1..n,1..m) general case
        for ( Sequence::const_iterator iX=bX; iX != eX; ++iX )
        {
            dScore = maxScorePrevRow + (*iX==baseY ? MATCH : MIS_MATCH);
            if ( pTrace[-1] == Trace::LEFT )
                hScore = std::max(maxScore+GAP_OPEN,hScore+GAP_EXTEND);
            else
                hScore = maxScore + GAP_OPEN;
            maxScorePrevRow = *pMaxScore;
            if ( pTrace[-nCols-1] == Trace::UP )
                vScore = std::max(maxScorePrevRow+GAP_OPEN,*pVScore+GAP_EXTEND);
            else
                vScore = maxScorePrevRow + GAP_OPEN;
            maxScore = max3(dScore,hScore,vScore,pTrace++);

            *pVScore++ = vScore;
            *pMaxScore++ = maxScore;
        }
    }
    return maxScore;
}

void SmithWatermanizer::traceback( Sequence::const_iterator bX,
                                   Sequence::const_iterator eX,
                                   Sequence::const_iterator bY,
                                   Sequence::const_iterator eY,
                                   std::ostream* pOS )
{
    using std::distance;
    auto nCols = distance(bX, eX);
    mMatchup.clear();
    if ( pOS )
    {
        mLineX.clear();
        mLineY.clear();
        mLineZ.clear();
    }
    Trace* pTB = &mTraceback.back();
    unsigned len = 0;
    while ( *pTB != Trace::DONE )
    {
        Trace trace = *pTB;
        if ( trace == Trace::UP )
        {
            if ( pOS )
            {
                mLineX.push_back(' ');
                mLineY.push_back(charFor(*--eY));
                mLineZ.push_back('|');
            }
            pTB -= nCols + 1;
        }
        else if ( trace == Trace::DIAG )
        {
            if ( pOS )
            {
                mLineX.push_back(charFor(*--eX));
                mLineY.push_back(charFor(*--eY));
                mLineZ.push_back(*eX==*eY?' ':'*');
            }
            pTB -= nCols + 2;
        }
        else
        {
            if ( pOS )
            {
                mLineX.push_back(charFor(*--eX));
                mLineY.push_back(' ');
                mLineZ.push_back('|');
            }
            pTB -= 1;
        }
        len += 1;
        if ( trace != *pTB )
        {
            mMatchup.emplace_back(trace,len);
            len = 0;
        }
    }
    std::reverse(mMatchup.begin(),mMatchup.end());
    if ( pOS )
    {
        if ( bX!=eX || bY!=eY )
            FatalErr("Internal error in traceback.");
        std::reverse(mLineX.begin(), mLineX.end());
        *pOS << mLineX;
        std::reverse(mLineY.begin(), mLineY.end());
        *pOS << mLineY;
        std::reverse(mLineZ.begin(), mLineZ.end());
        *pOS << mLineZ;
    }
}

// processes reads, looking for alignments that overlap the ORF, and counting up
// the codon values that we see.
class ORFCaller
{
public:
    ORFCaller( Reference const& ref, double minSW, Qual minQ,
                AminoAcidNamer const& aaNamer, std::ostream* pJunkOS )
    : mRef(ref), mMinPerBaseSWScore(minSW), mMinQ(minQ), mAANamer(aaNamer),
      mCCV(ref.getORFSize()/3), mNReads(0), mNAlignedReads(0),
      mNFPIndelReads(0), mNFSIndelReads(0), mpJunkOS(pJunkOS)
    {}

    void processRead( Sequence const& seq, Quals const& quals,
                        char const* fileName, size_t readId )
    { Alignment aln;
      mIndelLocs.clear();
      if ( findAlignment(seq,quals,fileName,readId,&aln) )
        processAlignment(aln,seq,quals);
      processIndelLocs(); }

    void processPairedRead( Sequence const& seq1, Quals const& quals1,
                            Sequence const& seq2, Quals const& quals2,
                            char const* fileName1, char const* fileName2,
                            size_t readId )
    { Alignment aln1, aln2;
      mIndelLocs.clear();
      if ( findAlignment(seq1,quals1,fileName1,readId,&aln1) )
      { if ( !findAlignment(seq2,quals2,fileName2,readId,&aln2) )
          processAlignment(aln1,seq1,quals1);
        else if ( Alignment::disjoint(aln1,aln2) )
        { processAlignment(aln1,seq1,quals1);
          processAlignment(aln2,seq2,quals2); }
        else if ( aln1.getRefStart() <= aln2.getRefStart() )
          processAlignmentPair(aln1,aln2,seq1,seq2,quals1,quals2);
        else
          processAlignmentPair(aln2,aln1,seq2,seq1,quals2,quals1); }
      else if ( findAlignment(seq2,quals2,fileName2,readId,&aln2) )
        processAlignment(aln2,seq2,quals2);
      processIndelLocs(); }

    double getFractionAligned() const
    { return mNReads ? 1.*mNAlignedReads/mNReads : 0.; }

    void writeCodonCounts( std::string const& filename );
    void writeAminoAcidCounts( std::string const& filename );
    void writeIndelCounts( std::string const& filename );
    void report( std::ostream& os );

private:
    bool findAlignment( Sequence const& seq, Quals const& quals,
                        char const* fileName, size_t readId, Alignment* pAln )
    { ++mNReads;
      if ( mRef.findLastAlignment(seq,quals,pAln) )
        return true;
      if ( !findIndels(*pAln,seq,quals) )
        dumpNonAligner(fileName,readId,seq,quals);
      return false; }

    bool findIndels( Alignment const& aln,
                        Sequence const& seq, Quals const& quals );
    void processIndelLocs();
    void processAlignment( Alignment const& aln,
                            Sequence const& seq, Quals const& quals );
    void processAlignmentPair( Alignment const& aln1, Alignment const& aln2,
                                Sequence const& seq1, Sequence const& seq2,
                                Quals const& quals1, Quals const& quals2 );
    void dumpNonAligner( char const* fileName, size_t readId,
                            Sequence const& seq, Quals const& quals )
    { if ( !mpJunkOS ) return;
      std::ostream& os = *mpJunkOS;
      os << '@' << fileName << ':' << readId << '\n';
      for ( Call call : seq ) os << charFor(call);
      os << "\n+\n";
      for ( Qual qual : quals ) os << (char)(qual+FASTQ_QUAL_OFFSET);
      os << '\n'; }

    struct IndelLoc
    {
        IndelLoc( Offset codonId, bool isFrameShift )
        : mCodonId(codonId), mFrameShift(isFrameShift) {}

        friend bool operator<( IndelLoc const& il1, IndelLoc const& il2 )
        { if ( il1.mCodonId != il2.mCodonId )
            return il1.mCodonId < il2.mCodonId;
          return il1.mFrameShift < il2.mFrameShift; }
        friend bool operator==( IndelLoc const& il1, IndelLoc const& il2 )
        { return il1.mCodonId == il2.mCodonId &&
                    il1.mFrameShift==il2.mFrameShift; }
        Offset mCodonId;
        bool mFrameShift;
    };

    Reference const& mRef;
    double mMinPerBaseSWScore;
    Qual mMinQ;
    AminoAcidNamer const& mAANamer;
    CodonCountsVec mCCV;
    size_t mNReads;
    size_t mNAlignedReads;
    size_t mNFPIndelReads;
    size_t mNFSIndelReads;
    std::ostream* mpJunkOS;
    Sequence mSeqTmp;
    Quals mQualsTmp;
    SmithWatermanizer mSW;
    std::vector<IndelLoc> mIndelLocs;
};


void ORFCaller::writeCodonCounts( std::string const& filename )
{
    std::ofstream ofs(filename);
    ofs <<
"AAA\tAAC\tAAG\tAAT\tACA\tACC\tACG\tACT\tAGA\tAGC\tAGG\tAGT\tATA\tATC\tATG\tATT\t"
"CAA\tCAC\tCAG\tCAT\tCCA\tCCC\tCCG\tCCT\tCGA\tCGC\tCGG\tCGT\tCTA\tCTC\tCTG\tCTT\t"
"GAA\tGAC\tGAG\tGAT\tGCA\tGCC\tGCG\tGCT\tGGA\tGGC\tGGG\tGGT\tGTA\tGTC\tGTG\tGTT\t"
"TAA\tTAC\tTAG\tTAT\tTCA\tTCC\tTCG\tTCT\tTGA\tTGC\tTGG\tTGT\tTTA\tTTC\tTTG\tTTT";

    for ( CodonCounts const& cCounts : mCCV )
    {
        char prefix = '\n';
        for ( size_t count : cCounts.mCounts )
        {
            ofs << prefix << count;
            prefix = '\t';
        }
    }
    ofs << '\n';
    ofs.close();
    if ( !ofs )
        FatalErr("Couldn't write codon counts to " << filename);
}

void ORFCaller::writeAminoAcidCounts( std::string const& filename )
{
    std::ofstream ofs(filename);
    char const* prefixS = "";
    for ( char label : mAANamer )
    {
        ofs << prefixS << label;
        prefixS = "\t";
    }

    std::vector<size_t> aaCounts;
    for ( CodonCounts const& cCounts : mCCV )
    {
        mAANamer.getAACounts(cCounts.mCounts,&aaCounts);
        char prefix = '\n';
        for ( size_t count : aaCounts )
        {
            ofs << prefix << count;
            prefix = '\t';
        }
    }
    ofs << '\n';
    ofs.close();
    if ( !ofs )
        FatalErr("Couldn't write amino acid counts to " << filename);
}

void ORFCaller::writeIndelCounts( std::string const& filename )
{
    std::ofstream ofs(filename);
    ofs << "NFS\tFS\tTotal\n";

    for ( CodonCounts const& cCounts : mCCV )
    {
        size_t total = cCounts.mIndelCounts[0]+cCounts.mIndelCounts[1];
        total = std::accumulate(cCounts.mCounts,cCounts.mCounts+NCODONS,total);
        ofs << cCounts.mIndelCounts[0] << '\t'
            << cCounts.mIndelCounts[1] << '\t'
            << total << '\n';
    }
    ofs.close();
    if ( !ofs )
        FatalErr("Couldn't write indel counts to " << filename);
}

void ORFCaller::report( std::ostream& os )
{
    if ( !mNAlignedReads )
        FatalErr("Not even a single read aligned.");

    Sequence orfSeq = mRef.getORFSeq();
    Sequence::const_iterator iORF = orfSeq.cbegin();
    unsigned codonId = 0;

    size_t minAA = std::numeric_limits<size_t>::max();
    size_t minAAIdx = 0;
    size_t maxAA = std::numeric_limits<size_t>::min();
    size_t maxAAIdx = 0;
    size_t totAA = 0;
    os << std::setprecision(5) << std::fixed;
    os << "\nCodon\tWT\tmin%sub\tmean%\tmax%\n";
    std::vector<size_t> aaCounts;
    for ( CodonCounts const& cCounts : mCCV )
    {
        mAANamer.getAACounts(cCounts.mCounts,&aaCounts);
        size_t aaId = mAANamer.indexFor(iORF);
        iORF += 3;
        size_t min = std::numeric_limits<size_t>::max();
        size_t max = std::numeric_limits<size_t>::min();
        size_t tot = 0;
        size_t aa = 0;
        auto iWildType=aaCounts.begin()+aaId;
        for ( auto itr=aaCounts.begin(),end=aaCounts.end(); itr != end; ++itr )
        {
            if ( itr == iWildType ) continue;
            size_t count = *itr;
            if ( count < min ) min = count;
            if ( count > max ) max = count;
            if ( count ) ++aa;
            tot += count;
        }
        double dTot = (tot+*iWildType)/100.;
        os << ++codonId
                  << '\t' << mAANamer[aaId]
                  << '\t' << min/dTot
                  << '\t' << tot/dTot/(aaCounts.size()-1)
                  << '\t' << max/dTot
                  << '\n';
        if ( aa < minAA ) { minAA = aa; minAAIdx = codonId; }
        if ( aa > maxAA ) { maxAA = aa; maxAAIdx = codonId; }
        totAA += aa;
    }
    os << "\nRepresented non-WT AA counts:\n";
    os << "min=" << minAA << " for codon " << minAAIdx << '\n';
    os << "mean=" << 1.*totAA/mCCV.size() << '\n';
    os << "max=" << maxAA << " for codon " << maxAAIdx << '\n';

    os << "\nIndel rates:\n";
    os << 100.*mNFPIndelReads/mNAlignedReads
         << "% of aligning reads contained frame-preserving indels.\n";
    os << 100.*mNFSIndelReads/mNAlignedReads
         << "% of aligning reads contained frame-shifting indels.\n";
}

bool ORFCaller::findIndels( Alignment const& aln,
                                Sequence const& seq, Quals const& )
{
    if ( aln.empty() )
        return false;

    Alignment alnBeg;
    if ( !mRef.findFirstAlignment(seq, &alnBeg)
            || aln.isReadRC() != alnBeg.isReadRC() )
        return false;

    Sequence const& ref = mRef.getRef();
    Sequence::const_iterator refBeg, refEnd, seqBeg, seqEnd;

    if ( !aln.isReadRC() )
    {
        if ( alnBeg.getRefStart() >= aln.getRefEnd() )
            return false;
        refBeg = ref.begin() + alnBeg.getRefStart();
        refEnd = ref.begin() + aln.getRefEnd();
        seqBeg = seq.begin() + alnBeg.getReadStart();
        seqEnd = seq.begin() + aln.getReadEnd();
    }
    else
    {
        mSeqTmp.clear();
        mSeqTmp.reserve(seq.size());
        for ( auto itr = seq.rbegin(), end = seq.rend(); itr != end; ++itr )
            mSeqTmp.push_back(complement(*itr));
        if ( aln.getRefStart() >= alnBeg.getRefEnd() )
            return false;
        refBeg = ref.begin() + aln.getRefStart();
        refEnd = ref.begin() + alnBeg.getRefEnd();
        seqBeg = mSeqTmp.begin() + aln.getReadStart();
        seqEnd = mSeqTmp.begin() + alnBeg.getReadEnd();
    }

    int score = mSW.align(seqBeg, seqEnd, refBeg, refEnd);
    double perBaseScore = 1. * score / (refEnd - refBeg);
    if ( perBaseScore < mMinPerBaseSWScore )
        return false;

    auto iExon = mRef.getORF().begin();
    using std::distance;
    bool hasFSIndel = false;
    bool hasNFSIndel = false;
    Offset refPos = distance(ref.begin(), refBeg);
    Offset orfPos = 0;
    for ( SmithWatermanizer::Block const& block : mSW.getMatchup() )
    {
        while ( iExon->mEnd <= refPos )
        {
            orfPos += iExon->size();
            ++iExon;
        }
        bool fs = block.mLen % 3;
        if ( block.mDir != SmithWatermanizer::Trace::DIAG
                && refPos >= iExon->mBegin )
        {
            mIndelLocs.emplace_back((orfPos+refPos-iExon->mBegin)/3,fs);
            if ( fs ) hasFSIndel = true;
            else hasNFSIndel = true;
        }
        refPos += block.mLen;
    }
    if ( hasFSIndel )
        mNFSIndelReads += 1;
    else if ( hasNFSIndel )
        mNFPIndelReads += 1;
    ++mNAlignedReads;
    return true;
}

void ORFCaller::processIndelLocs()
{
    if ( mIndelLocs.size() )
    {
        auto beg = mIndelLocs.begin(), end = mIndelLocs.end();
        std::sort(beg,end);
        mIndelLocs.erase(std::unique(beg,end),end);
        for ( IndelLoc const& il : mIndelLocs )
            mCCV[il.mCodonId].mIndelCounts[il.mFrameShift] += 1;
    }
}

void ORFCaller::processAlignment( Alignment const& alnA,
                                    Sequence const& seq, Quals const& quals )
{
    ++mNAlignedReads;

    Alignment aln(alnA);
    auto oItr = mRef.getORF().begin();
    Offset orfStart = 0;
    while ( oItr->mEnd <= aln.getRefStart() )
    {
        orfStart += oItr->size();
        ++oItr;
    }

    if ( aln.getRefStart() > oItr->mBegin )
        orfStart += aln.getRefStart() - oItr->mBegin;
    CodonAccumulator acc(mCCV.begin(),orfStart,mMinQ);

    while ( oItr->mBegin < aln.getRefEnd() )
    {
        if ( oItr->mBegin > aln.getRefStart() )
            aln += oItr->mBegin - aln.getRefStart();

        Offset len = std::min(aln.getRefEnd(),oItr->mEnd)-aln.getRefStart();
        if ( aln.isReadRC() )
        {
            auto qItr = quals.rbegin()+aln.getReadStart();
            auto sItr = seq.rbegin()+aln.getReadStart();
            for ( auto sEnd = sItr+len; sItr != sEnd; ++sItr,++qItr )
                acc(complement(*sItr),*qItr);
        }
        else
        {
            auto qItr = quals.begin()+aln.getReadStart();
            auto sItr = seq.begin()+aln.getReadStart();
            for ( auto sEnd = sItr+len; sItr != sEnd; ++sItr,++qItr )
                acc(*sItr,*qItr);
        }
        aln += len;
        if ( aln.empty() )
            break;
        ++oItr;
    }
}

// make a super-read out of the pair.
// this is super ugly.  sorry.  can't figure out anything clean and elegant.
void ORFCaller::processAlignmentPair( Alignment const& aln1,
        Alignment const& aln2, Sequence const& seq1, Sequence const& seq2,
        Quals const& quals1, Quals const& quals2 )
{
    if ( aln1.isReadRC() == aln2.isReadRC() )
        return; // one of the alignments is bogus
                // pairs cannot be on the same strand

    if ( aln1.getRefStart() > aln2.getRefStart() )
        FatalErr("Internal error: aln1 starts after aln2.");

    mSeqTmp.clear();
    mSeqTmp.reserve(seq1.size()+seq2.size());
    mQualsTmp.clear();
    mQualsTmp.reserve(quals1.size()+quals2.size());

    if ( aln1.isReadRC() )
    {
        for ( auto itr=seq1.rbegin()+aln1.getReadStart(),
                end=seq1.rbegin()+aln1.getReadEnd(); itr != end; ++itr )
            mSeqTmp.push_back(complement(*itr));
        mQualsTmp.assign(quals1.rbegin()+aln1.getReadStart(),
                            quals1.rbegin()+aln1.getReadEnd());
        Offset overlapLen = std::min(aln1.getRefEnd(),aln2.getRefEnd())-
                                aln2.getRefStart();
        Offset start1 = aln2.getRefStart() - aln1.getRefStart();
        auto s1 = mSeqTmp.begin()+start1;
        auto s2 = seq2.begin()+aln2.getReadStart();
        auto q2 = quals2.begin()+aln2.getReadStart();
        for ( auto q1=mQualsTmp.begin()+start1,
                   qEnd=q1+overlapLen; q1 != qEnd; ++s1,++q1,++s2,++q2 )
            if ( *q2 > *q1 ) { *s1 = *s2; *q1 = *q2; }
        if ( aln2.getRefEnd() > aln1.getRefEnd() )
        {
            Offset extra = aln2.getRefEnd()-aln1.getRefEnd();
            auto sEnd = seq2.begin()+aln2.getReadEnd();
            std::copy(sEnd-extra,sEnd,std::back_inserter(mSeqTmp));
            auto qEnd = quals2.begin()+aln2.getReadEnd();
            std::copy(qEnd-extra,qEnd,std::back_inserter(mQualsTmp));
        }
    }
    else
    {
        mSeqTmp.assign(seq1.begin()+aln1.getReadStart(),
                        seq1.begin()+aln1.getReadEnd());
        mQualsTmp.assign(quals1.begin()+aln1.getReadStart(),
                            quals1.begin()+aln1.getReadEnd());
        Offset overlapLen = std::min(aln1.getRefEnd(),aln2.getRefEnd())-
                                aln2.getRefStart();
        Offset start1 = aln2.getRefStart() - aln1.getRefStart();
        auto s1 = mSeqTmp.begin()+start1;
        auto s2 = seq2.rbegin()+aln2.getReadStart();
        auto q2 = quals2.rbegin()+aln2.getReadStart();
        for ( auto q1=mQualsTmp.begin()+start1,
                   qEnd=q1+overlapLen; q1 != qEnd; ++s1,++q1,++s2,++q2 )
            if ( *q2 > *q1 ) { *s1 = complement(*s2); *q1 = *q2; }
        if ( aln2.getRefEnd() > aln1.getRefEnd() )
        {
            Offset extra = aln2.getRefEnd()-aln1.getRefEnd();
            auto sEnd = seq2.rbegin()+aln2.getReadEnd();
            for ( auto sItr = sEnd-extra; sItr != sEnd; ++sItr )
                mSeqTmp.push_back(complement(*sItr));
            auto qEnd = quals2.rbegin()+aln2.getReadEnd();
            std::copy(qEnd-extra,qEnd,std::back_inserter(mQualsTmp));
        }
    }
    Alignment aln(aln1.getRefStart(),
                    std::max(aln1.getRefEnd(),aln2.getRefEnd()),
                    0,false);
    processAlignment( aln, mSeqTmp, mQualsTmp );
    ++mNAlignedReads;
}



void processFastqs( char** beg, char** end, ORFCaller* pCaller )
{
    std::vector<Sequence> seqs;
    std::vector<Quals> quals;
    size_t fileId = 0;
    while ( beg != end )
    {
        char const* fileName = *beg++;
        std::cerr << "Reading " << fileName << "..." << std::endl;
        FastqReader fqr(fileName);
        size_t readId = 0;
        while ( fqr.getChunk(seqs,quals) )
        {
            auto qItr = quals.begin();
            for ( Sequence const& sq : seqs )
            {
                pCaller->processRead(sq,*qItr,fileName,readId++);
                ++qItr;
            }
        }
        std::cerr << readId << " reads in file " << fileId << ".\n";
        fileId += 1;
    }
}

void processPairedFastqs( char** beg, char** end, ORFCaller* pCaller )
{
    if ( (beg-end) & 1 )
        FatalErr("You requested paired mode, but there are an "
                 "odd number of fastq files.");

    std::vector<Sequence> seqs1, seqs2;
    std::vector<Quals> quals1, quals2;
    size_t fileId = 0;
    while ( beg != end )
    {
        char const* fileName1 = *beg++;
        char const* fileName2 = *beg++;
        std::cerr << "Reading " << fileName1 << " and " << fileName2 << "..."
                    << std::endl;
        FastqReader fqr1(fileName1);
        FastqReader fqr2(fileName2);

        size_t readId = 0;
        while ( fqr1.getChunk(seqs1,quals1) )
        {
            size_t nReads = seqs1.size();
            if ( !fqr2.getChunk(seqs2,quals2) || nReads != seqs2.size() )
                FatalErr("Unequal numbers of reads in the paired fastq files.");

            auto sItr1=seqs1.begin(), sEnd=seqs1.end(), sItr2=seqs2.begin();
            auto qItr1=quals1.begin(), qItr2=quals2.begin();
            for ( ; sItr1 != sEnd; ++sItr1,++sItr2,++qItr1,++qItr2 )
                pCaller->processPairedRead(*sItr1,*qItr1,*sItr2,*qItr2,
                                            fileName1,fileName2,readId++);
        }
        if ( fqr2.getChunk(seqs2,quals2) )
            FatalErr("Unequal numbers of reads in the paired fastq files.");
        std::cerr << readId << " reads in file pair " << fileId << ".\n";
        fileId += 1;
    }
}

} // end of anonymous namespace

int main( int argc, char** argv )
{
    bool dumpNonAligners = false;
    unsigned K = 25;
    double maxMismatchFrac = .1;
    std::string outputFilename = "orfcall";
    bool pairedMode = false;
    unsigned minQ = 20;
    double minPerBaseSWScore = 2.;
    std::string codonTranslation =
            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF";

    int opt;
    while ( (opt = getopt(argc, argv, "dK:M:O:pQ:S:T:")) != -1 )
    {
        switch ( opt )
        {
        case 'd':
            dumpNonAligners = true;
            break;
        case 'K':
          { std::istringstream iss(optarg);
            if ( !(iss >> K) || K < 14 || K > 64 )
                FatalErr("Can't interpret K='" << optarg
                        << "'. Expecting an integer between 14 and 64"
                           " for kmer size."); }
            // note to self:  we don't have to worry about palindromes because
            // we're adding the rc of the ref to the dictionary, too, and that
            // will cause palindromes to be discarded as dups.
            break;
        case 'M':
          { std::istringstream iss(optarg);
            if ( !(iss >> maxMismatchFrac) ||
                    maxMismatchFrac < 0. || maxMismatchFrac > 1. )
                FatalErr("Can't interpret M='" << optarg
                        << "'. Expecting a fraction between "
                           "0. and 1. for max mismatch rate."); }
            break;
        case 'O':
            outputFilename = optarg;
            break;
        case 'p':
            pairedMode = true;
            break;
        case 'Q':
          { std::istringstream iss(optarg);
            if ( !(iss >> minQ) || minQ > 60 )
                FatalErr("Can't interpret Q='" << optarg
                        << "'. Expecting an integer min quality score."); }
            break;
        case 'S':
          { std::istringstream iss(optarg);
            if ( !(iss >> minPerBaseSWScore) )
                FatalErr("Can't interpret S='" << optarg
                        << "'. Expecting a minimum per-base S/W score."); }
            break;
        case 'T':
            codonTranslation = optarg;
            if ( codonTranslation.size() != NCODONS )
                FatalErr("Can't interpret T='" << optarg <<
                        " as a table of codon calls.  Expected 64 characters.");
            break;
        case '?':
            if ( optopt == 'K' )
                FatalErr("Need value for K (kmer size) option.");
            else if ( optopt == 'M' )
                FatalErr("Need value for M (max mismatch fraction) option.");
            else if ( optopt == 'O' )
                FatalErr("Need value for O (output filename) option.");
            else if ( optopt == 'Q' )
                FatalErr("Need value for Q (min quality score) option.");
            else if ( optopt == 'S' )
                FatalErr("Need value for S (min per-base S/W score.");
            else if ( optopt == 'T' )
                FatalErr("Need value for T (table of codon calls) option.");
            else
                FatalErr("Unknown option character " << optopt);
            break;
        }
    }

    if ( optind >= argc )
        FatalErr("Usage: ORFCall [-d] [-K KMER_SIZE] [-M MAX_MISMATCH_FRACTION] "
                 "[-O OUTPUT_FILENAME] [-p] [-Q MIN_QUAL] [-T CODON_TABLE] "
                 "ref.fasta reads1.fastq [reads2.fastq...]\n"
                 "The -p flag signals paired mode.\n");

    Reference ref(argv[optind++],K,minQ,maxMismatchFrac);
    ref.report();

    std::ofstream* pJunkOS = nullptr;
    if ( dumpNonAligners )
        pJunkOS = new std::ofstream(outputFilename+".unaligned.fasta");
    AminoAcidNamer aaNamer(codonTranslation);
    ORFCaller caller(ref,minPerBaseSWScore,minQ,aaNamer,pJunkOS);
    char** beg = argv+optind;
    char** end = argv+argc;
    if ( !pairedMode )
        processFastqs(beg,end,&caller);
    else
        processPairedFastqs(beg,end,&caller);
    std::cerr << "Aligned reads: " << 100.*caller.getFractionAligned() << "%\n";

    if ( pJunkOS )
        delete pJunkOS;
    caller.writeCodonCounts(outputFilename+".cdn");
    caller.writeAminoAcidCounts(outputFilename+".aa");
    caller.writeIndelCounts(outputFilename+".indels");
    caller.report(std::cout);
}
