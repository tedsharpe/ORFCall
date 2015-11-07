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

// a 6-bit number designating one of the 64 possible codons
// (AAA=0, AAC=1, ... TTT=63)
typedef unsigned Codon;
Codon const NCODONS = 64;
std::string seqFor( Codon codon )
{
    std::string result;
    result.reserve(3);
    result.push_back("ACGT"[(codon >> 4) & 3]);
    result.push_back("ACGT"[(codon >> 2) & 3]);
    result.push_back("ACGT"[codon & 3]);
    return result;
}

typedef unsigned char Call;
Call const NCALL = 4;
inline Call complement( Call call ) { return (call == 4) ? 4 : (call ^ 3); }
inline char charFor( Call call ) { return "ACGTN"[call]; }

typedef std::vector<Call> Sequence;
typedef Sequence::size_type Offset;
typedef unsigned char Qual;
typedef std::vector<Qual> Quals;

// describes an exon as a pair of offsets on the reference sequence.
struct Exon
{
    Offset size() const { return mEnd - mBegin; }

    Offset mBegin;
    Offset mEnd;
};
typedef std::vector<Exon> ORF;

typedef Sequence::const_iterator Kmer;
inline Codon codonFor( Kmer kmer ) { return (kmer[0]<<4)|(kmer[1]<<2)|kmer[2]; }

// produces a hash value for the next K bases from a sequence.
class KmerHasher
{
public:
    KmerHasher( unsigned K ) : mK(K) {}

    size_t operator()( Kmer itr ) const
    {
        size_t result = 14695981039346656037ul;
        for ( Kmer end(itr + mK); itr != end; ++itr )
            result = 1099511628211ul * (result ^ *itr); // FNV-1a algorithm
        return result;
    }

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
    enum class Status : char { UNSET, UNIQUE, DUP };

public:
    KmerLoc() : mOffset(0), mIsRC(false), mStatus(Status::UNSET) {}

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
    Reference( char const* refFile, char const* planFile, unsigned const K )
    : mK(K), mUniqueKmers(0,KmerHasher(K),KmerComparator(K))
    {
        readRefSeq(refFile);
        readPlan(planFile);
        buildDictionary();
    }

    unsigned getK() const { return mK; }
    Sequence const& getRef() const { return mRefSeq; }

    ORF const& getORF() const { return mORF; }
    std::vector<Codon> const& getORFWTCodons() const { return mWTCodons; }
    std::vector<int64_t> const& getPlannedMasks() const { return mPlannedMasks; }

    // finds a maximum-length correspondence between the reference sequence
    // and another sequence (or its RC) as evidenced by a matching kmer.
    // if no matching kmer is found, the Alignment will be empty.
    Alignment findAlignment( Sequence const& seq ) const
    {
        if ( seq.size() < mK )
            return Alignment();

        auto notFound = mUniqueKmers.end();
        for ( Kmer kmer(seq.begin()), end(seq.end()-mK+1); kmer != end; ++kmer )
        {
            auto itr = mUniqueKmers.find(kmer);
            if ( itr != notFound )
            {
                Offset refStart = itr->second.getOffset();
                Offset readStart = kmer - seq.begin();
                Offset backOff = std::min(refStart,readStart);
                refStart -= backOff;
                readStart -= backOff;
                Offset len = std::min(mRefSeq.size()-refStart,
                                        seq.size()-readStart);
                bool isRC = itr->second.isRC();
                if ( isRC )
                {
                    refStart = mRefSeq.size()-refStart-len;
                    readStart = seq.size()-readStart-len;
                }
                return Alignment(refStart,refStart+len,readStart,isRC);
            }
        }
        return Alignment();
    }

    void report() const
    {
        size_t refLen = mRefSeq.size();
        std::cerr << "Reference length: " << refLen << ".\n";
        std::cerr << "Dict size: " << mUniqueKmers.size() << " ("
                << 2 * (refLen - mK + 1) - mUniqueKmers.size()
                << " duplicate kmers).\n";
        std::cerr << "ORF length in codons: " << mWTCodons.size() << ".\n";
        std::cerr << "ORF segments:";
        char prefix = ' ';
        // leave final sentinel out of description
        for ( auto itr = mORF.begin(), end = mORF.end()-1; itr != end; ++itr )
        {
            std::cerr << prefix << itr->mBegin + 1 << ':' << itr->mEnd + 1;
            prefix = ',';
        }
        std::cerr << ".\n";
    }

private:
    void readRefSeq( char const* refFile );
    void readPlan( char const* planFile );
    void buildDictionary();

    unsigned mK;
    Sequence mRefSeq;
    Sequence mRefSeqRC;
    KmerDict mUniqueKmers;

    // ref seq boundaries for the ORF (i.e., exon locations on ref)
    ORF mORF;

    // WT codon value for each codon in the ORF
    std::vector<Codon> mWTCodons;

    // bit mask of other codon values that we might expect to see for each codon
    std::vector<int64_t> mPlannedMasks;
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
                inORF=true; mORF.push_back(Exon{mRefSeq.size(),0}); break;
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
        mORF.push_back(Exon{0,mRefSeq.size()});
    }
    if ( !ifs.eof() )
        FatalErr("Failed to read reference file to the end.");

    if ( mRefSeq.size() < mK )
        FatalErr("There are only " << mRefSeq.size()
                    << "bases in the reference file.  K is too big.");

    // add sentinel
    mORF.push_back(Exon{mRefSeq.size(),mRefSeq.size()});

    size_t orfSize = 0ul;
    for ( Exon const& exon : mORF )
        orfSize += exon.size();
    if ( orfSize % 3 )
        FatalErr("There are " << orfSize <<
                " bases in the ORF, which isn't divisible by 3.");

    Sequence orfSeq;
    orfSeq.reserve(orfSize);
    for ( Exon const& exon : mORF )
        std::copy(mRefSeq.begin()+exon.mBegin,mRefSeq.begin()+exon.mEnd,
                    std::back_inserter(orfSeq) );

    size_t nCodons = orfSize / 3;
    mWTCodons.reserve(nCodons);
    for ( Kmer itr(orfSeq.begin()), end(orfSeq.end()); itr < end; itr += 3 )
        mWTCodons.push_back(codonFor(itr));
}

void Reference::readPlan( char const* planFile ) {
    std::ifstream ifs(planFile,std::ios_base::in|std::ios_base::binary);
    if ( !ifs )
        FatalErr("Can't open plan file: " << planFile);

    std::string line;
    if ( !std::getline(ifs,line) || line.size() != 64*3+63 )
        FatalErr("Can't read plan file header.");

    size_t nCodons = mWTCodons.size();
    size_t lineNo = 0;
    mPlannedMasks.reserve(nCodons);
    while ( std::getline(ifs,line) )
    {
        ++lineNo;
        if ( line.size() != 64+63 )
            FatalErr("Plan file line " << lineNo << " is weird:\n" << line);
        int64_t val = 0;
        for ( Codon codon = 0; codon < NCODONS-1; ++codon )
        {
            char chr = line[2*codon+1];
            if ( chr != '\t' && chr != ' ' )
                FatalErr("Plan file line " << lineNo << " is weird:\n" << line);
        }
        for ( Codon codon = 0; codon < NCODONS; ++codon )
        {
            switch (line[2*codon])
            {
            case '0': break;
            case '1': val |= 1ul << codon; break;
            case '2':
                if ( lineNo <= nCodons )
                {
                    Codon refCodon = mWTCodons[lineNo-1];
                    if ( refCodon != codon )
                        FatalErr("Plan file line " << lineNo <<
                                " indicates WT for codon " << seqFor(codon) <<
                                " but the reference sequence has codon "
                                << seqFor(refCodon) << " at that position.");
                }
                break;
            default:
                FatalErr("Plan file line " << lineNo << " is weird:\n" << line);
            }
        }
        mPlannedMasks.push_back(val);
    }
    if ( mPlannedMasks.size() > nCodons )
        FatalErr("Plan file has " << mPlannedMasks.size() <<
                    " lines, but there are only " << nCodons << " codons.");
    else if ( mPlannedMasks.size() < nCodons )
    {
        std::cerr << "Warning: Plan file has only " << mPlannedMasks.size() <<
            " lines, but there are " << nCodons << " codons." << std::endl;
        mPlannedMasks.resize(nCodons);
    }
}

void Reference::buildDictionary()
{
    mRefSeqRC.reserve(mRefSeq.size());
    for ( auto itr=mRefSeq.rbegin(),end=mRefSeq.rend(); itr != end; ++itr )
        mRefSeqRC.push_back(complement(*itr));

    // build a dictionary of unique kmers in ref and ref RC
    mUniqueKmers.reserve(2*mRefSeq.size());
    Offset offset = 0;
    Kmer end(mRefSeq.end()-mK+1);
    for ( Kmer itr(mRefSeq.begin()); itr != end; ++itr )
        mUniqueKmers[itr].set(offset++,false);

    offset = 0;
    end = mRefSeqRC.end()-mK+1;
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
    {
        auto len = strlen(fastqFile);
        auto flags = std::ios_base::in | std::ios_base::binary;
        if ( len >= 3 && !memcmp(fastqFile + len - 3, ".gz", 3) )
            return new igzstream(fastqFile, flags);
        return new std::ifstream(fastqFile, flags);
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

// a translator of codons (6-bit values encoding 3 base calls) into amino acid symbols.
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
      mLabels.resize(size_t(std::unique(beg,end)-beg));
      end = mLabels.end();
      for ( size_t idx = 0; idx != NCODONS; ++idx )
        mCodonToAA[idx] = std::lower_bound(beg,end,codonTranslation[idx])-beg; }

    void getAACounts( size_t const* codonCounts,
                            std::vector<size_t>* pAACounts ) const
    { std::vector<size_t>& aaCounts = *pAACounts;
      aaCounts.clear(); aaCounts.resize(mLabels.size());
      for ( size_t idx : mCodonToAA )
          aaCounts[idx] += *codonCounts++; }

    // translate a 6-bit codon value to an amino acid index
    size_t indexFor( Codon codon ) const { return mCodonToAA[codon]; }

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
    CodonCounts() : mCounts{}, mIndelCounts{}, mCoverage(0) {}

    size_t mCounts[NCODONS];
    size_t mIndelCounts[2];
    size_t mCoverage;
};
typedef std::vector<CodonCounts> CodonCountsVec;
typedef CodonCountsVec::iterator CCVItr;

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

    std::vector<Block> const& align( Sequence::const_iterator bX,
                                     Sequence::const_iterator eX,
                                     Sequence::const_iterator bY,
                                     Sequence::const_iterator eY,
                                     std::ostream* pOS=nullptr )
    { score(bX,eX,bY,eY);
      traceback(bX,eX,bY,eY,pOS);
      return mMatchup; }

private:
    static int const GAP_OPEN = -14;
    static int const GAP_EXTEND = -3;
    static int const MATCH = 4;
    static int const MIS_MATCH = -6;

    void score( Sequence::const_iterator bX, Sequence::const_iterator eX,
                Sequence::const_iterator bY, Sequence::const_iterator eY );

    void traceback( Sequence::const_iterator bX, Sequence::const_iterator eX,
                    Sequence::const_iterator bY, Sequence::const_iterator eY,
                    std::ostream* pOS=nullptr );

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

    // where to start the traceback
    Trace* mpBestTracebackStart;

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

void SmithWatermanizer::score( Sequence::const_iterator bX,
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
    int rowInit = 0;//GAP_OPEN;
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

    int bestRowScore = maxScore;
    Trace* pBestTrace = pTrace;

    int maxScorePrevRow;
    while ( ++bY != eY )
    {
        baseY = *bY;
        pVScore = &mVScores.front();
        pMaxScore = &mMaxScores.front();
        *pTrace++ = Trace::UP;

        // (0,1..m) can't extend hGap
        dScore = rowInit + (baseX0==baseY ? MATCH : MIS_MATCH);
        //rowInit += GAP_EXTEND;
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

        if ( maxScore > bestRowScore )
        {
            bestRowScore = maxScore;
            pBestTrace = pTrace;
        }
    }
    mpBestTracebackStart = pBestTrace-1;
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
    Trace* pTB = mpBestTracebackStart;
    eY -= (&mTraceback.back() - pTB - 1)/(nCols + 1);

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

        Block const& firstBlock = mMatchup[0];
        if ( firstBlock.mDir == Trace::UP )
        {
            len = mLineX.size() - firstBlock.mLen;
            mLineX.resize(len);
            mLineY.resize(len);
            mLineZ.resize(len);
        }
        std::reverse(mLineX.begin(), mLineX.end());
        *pOS << mLineX << '\n';
        std::reverse(mLineY.begin(), mLineY.end());
        *pOS << mLineY << '\n';
        std::reverse(mLineZ.begin(), mLineZ.end());
        *pOS << mLineZ << '\n';
    }
}

enum class VariationType { ALTERNATE_CODON, INSERTION, DELETION };
struct Variation
{
    VariationType mType;
    Offset mCodonIdx;
    // when mType is ALTERNATE_CODON, mValue is the alternate Codon
    // when INSERTION or DELETION, mValue is the length of the indel
    unsigned mValue;

    friend bool operator<( Variation const& var1, Variation const& var2 )
    { if ( var1.mType != var2.mType ) return var1.mType < var2.mType;
      if ( var1.mCodonIdx != var2.mCodonIdx )
          return var1.mCodonIdx < var2.mCodonIdx;
      return var1.mValue < var2.mValue; }

    friend bool operator==( Variation const& var1, Variation const& var2 )
    { return var1.mType == var2.mType &&
            var1.mCodonIdx == var2.mCodonIdx &&
            var1.mValue == var2.mValue; }
};

struct ReadReport
{
    Offset mBeginCoverage;
    Offset mEndCoverage;
    size_t mNMismatches; // exclusive of those in the 1st alternate codon
    std::vector<Variation> mVariations;
};

class AlignmentScorer
{
public:
    AlignmentScorer( Reference const& ref, Offset refStart, Qual minQ )
    : mRefPos(refStart), mExonItr(ref.getORF().begin()), mNextExon{0,0},
      mWTBeg(ref.getORFWTCodons().begin()), mMinQ(minQ),
      mCodon(0), mPhase(0), mNCodonMismatches(0), mNHiQCalls(0),
      mFirstCodonCovered(0), mNMismatches(0)
    {
        while ( mExonItr->mEnd <= mRefPos )
        {
            mFirstCodonCovered += mExonItr->size();
            ++mExonItr;
        }
        mNextExon = *mExonItr;
        if ( mRefPos > mExonItr->mBegin )
            mFirstCodonCovered += mRefPos - mExonItr->mBegin;
        mPhase = mFirstCodonCovered % 3;
        mFirstCodonCovered = mFirstCodonCovered / 3;
        mWTItr = ref.getORFWTCodons().begin()+mFirstCodonCovered;
        mPlannedItr = ref.getPlannedMasks().begin()+mFirstCodonCovered;
        mCodon = *mWTItr >> ((3 - mPhase)*2);
    }

    void operator()( Call wtCall, Call call, Qual qual )
    {
        if ( mRefPos == mNextExon.mEnd )
            mNextExon = *++mExonItr;
        if ( mRefPos++ < mNextExon.mBegin )
        {
            if ( wtCall != call && qual >= mMinQ )
                mNMismatches += 1;
        }
        else
        {
            if ( qual >= mMinQ )
            {
                mNHiQCalls += 1;
                if ( wtCall != call )
                {
                    mNMismatches += 1;
                    mNCodonMismatches += 1;
                }
            }
            mCodon = (mCodon << 2) | call;
            if ( ++mPhase == 3 )
            {
                if ( mCodon != *mWTItr &&
                     ((1ul << mCodon) & *mPlannedItr) &&
                     mNHiQCalls == 3 )
                {
                    if ( mVariations.empty() )
                        mNMismatches -= mNCodonMismatches;
                    mVariations.push_back(
                        Variation{VariationType::ALTERNATE_CODON,
                                   Offset(mWTItr-mWTBeg),
                                   mCodon} );
                }
                mCodon = 0;
                mPhase = 0;
                mNCodonMismatches = 0;
                mNHiQCalls = 0;
                ++mWTItr;
                ++mPlannedItr;
            }
        }
    }

    bool inORF() const { return mRefPos >= mNextExon.mBegin; }
    Offset getCodonIdx() const { return Offset(mWTItr-mWTBeg); }

    ReadReport getReport() const
    {
        return ReadReport{mFirstCodonCovered,
                          Offset(mWTItr-mWTBeg),
                          mNMismatches,
                          mVariations};
    }

private:
    // stuff to track where we are in the ORF
    Offset mRefPos;
    ORF::const_iterator mExonItr;
    Exon mNextExon;
    std::vector<Codon>::const_iterator mWTBeg;
    std::vector<Codon>::const_iterator mWTItr;
    std::vector<int64_t>::const_iterator mPlannedItr;

    // a quality score must be at least this high to count a mismatch
    Qual mMinQ;

    // stuff to track the current codon
    Codon mCodon;
    unsigned mPhase;
    unsigned mNCodonMismatches;
    unsigned mNHiQCalls;

    // stuff to report at the end
    Offset mFirstCodonCovered;
    size_t mNMismatches;
    std::vector<Variation> mVariations;
};

// processes reads, looking for alignments that overlap the ORF, and counting up
// reads that contain a single alternate codon (according to the plan) or that
// contain indels.
class ORFCaller
{
public:
    ORFCaller( Reference const& ref, double minSW, Qual minQ,
                        unsigned maxSeqErrs, AminoAcidNamer const& aaNamer,
                        std::ostream* pJunkOS )
    : mRef(ref), mMinPerBaseSWScore(minSW), mMinQ(minQ),
      mMaxSeqErrs(maxSeqErrs), mAANamer(aaNamer),
      mCCV(ref.getORFWTCodons().size()), mNReads(0), mNAlignedReads(0),
      mNFPIndelReads(0), mNFSIndelReads(0), mpJunkOS(pJunkOS)
    {}

    void processFastqs( char** beg, char** end );

    void processPairedFastqs( char** beg, char** end );

    double getFractionAligned() const
    { return mNReads ? 1.*mNAlignedReads/mNReads : 0.; }

    void writeCodonCounts( std::string const& filename );
    void writeAminoAcidCounts( std::string const& filename );
    void writeIndelCounts( std::string const& filename );
    void report( std::ostream& os );

private:
    void processRead( Sequence const& seq, Quals const& quals,
                        char const* fileName, size_t readId )
    {
        mNReads += 1;
        Alignment aln = mRef.findAlignment(seq);
        if ( aln.empty() )
            dumpNonAligner(fileName,readId,seq,quals);
        else
            processAlignment(seq,quals,aln);
    }

    void processPairedRead( Sequence const& seq1, Quals const& quals1,
                            Sequence const& seq2, Quals const& quals2,
                            char const* fileName1, char const* fileName2,
                            size_t readId )
    {
        mNReads += 2;
        Alignment aln1 = mRef.findAlignment(seq1);
        Alignment aln2 = mRef.findAlignment(seq2);
        if ( aln1.empty() )
        {
            dumpNonAligner(fileName1,readId,seq1,quals1);
            if ( aln2.empty() )
                dumpNonAligner(fileName2,readId,seq2,quals2);
            else
                processAlignment(seq2,quals2,aln2);
        }
        else if ( aln2.empty() )
        {
            dumpNonAligner(fileName2,readId,seq2,quals2);
            processAlignment(seq1,quals1,aln1);
        }
        else if ( Alignment::disjoint(aln1,aln2) )
        {
            processAlignment(seq1,quals1,aln1);
            processAlignment(seq2,quals2,aln2);
        }
        else
        {
            mNAlignedReads += 2;
            ReadReport report1 = findAlternateCodons(seq1,quals1,aln1);
            if ( report1.mNMismatches > mMaxSeqErrs ||
                    report1.mVariations.size() > 1 )
                report1 = findIndels(aln1,seq1,quals1);
            ReadReport report2 = findAlternateCodons(seq2,quals2,aln2);
            if ( report2.mNMismatches > mMaxSeqErrs ||
                    report2.mVariations.size() > 1 )
                report2 = findIndels(aln2,seq2,quals2);
            if ( report1.mNMismatches <= mMaxSeqErrs )
            {
                if ( report2.mNMismatches > mMaxSeqErrs )
                    processReport(report1);
                else
                {
                    ReadReport combo{
                        std::min(report1.mBeginCoverage,report2.mBeginCoverage),
                        std::max(report1.mEndCoverage,report2.mEndCoverage),
                        report1.mNMismatches+report2.mNMismatches,
                        report1.mVariations
                    };
                    combo.mVariations.insert(combo.mVariations.end(),
                                                report2.mVariations.begin(),
                                                report2.mVariations.end());
                    auto beg(combo.mVariations.begin());
                    auto end(combo.mVariations.end());
                    std::sort(beg,end);
                    combo.mVariations.erase(std::unique(beg,end),end);
                    processReport(combo);
                }
            }
            else if ( report2.mNMismatches <= mMaxSeqErrs )
                processReport(report2);
        }
    }

    void processAlignment( Sequence const& seq, Quals const& quals,
                                Alignment const& aln )
    {
        mNAlignedReads += 1;
        ReadReport readReport = findAlternateCodons(seq,quals,aln);
        if ( readReport.mNMismatches > mMaxSeqErrs ||
                readReport.mVariations.size() > 1 )
            readReport = findIndels(aln,seq,quals);
        if ( readReport.mNMismatches <= mMaxSeqErrs )
            processReport(readReport);
    }

    ReadReport findAlternateCodons( Sequence const& seq,
                                    Quals const& quals,
                                    Alignment const& aln )
    {
        AlignmentScorer alignmentScorer(mRef,aln.getRefStart(),mMinQ);
        auto refItr(mRef.getRef().begin()+aln.getRefStart());
        auto refEnd(mRef.getRef().begin()+aln.getRefEnd());
        if ( !aln.isReadRC() )
        {
            auto readItr(seq.begin()+aln.getReadStart());
            auto qualItr(quals.begin()+aln.getReadStart());
            while ( refItr != refEnd )
            {
                alignmentScorer(*refItr,*readItr,*qualItr);
                ++refItr; ++readItr; ++qualItr;
            }
        }
        else
        {
            auto readItr(seq.rbegin()+aln.getReadStart());
            auto qualItr(quals.rbegin()+aln.getReadStart());
            while ( refItr != refEnd )
            {
                alignmentScorer(*refItr,complement(*readItr),*qualItr);
                ++refItr; ++readItr; ++qualItr;
            }
        }
        return alignmentScorer.getReport();
    }

    ReadReport findIndels( Alignment const& aln,
                            Sequence const& seq,
                            Quals const& quals );

    void processReport( ReadReport const& readReport )
    {
        for ( auto itr(mCCV.begin()+readReport.mBeginCoverage),
                end(mCCV.begin()+readReport.mEndCoverage);
                itr != end; ++itr )
            itr->mCoverage += 1;

        bool readHasFrameShiftIndel = false;
        bool readHasNonFrameShiftIndel = false;
        for ( Variation const& var : readReport.mVariations )
        {
            switch ( var.mType )
            {
            case VariationType::ALTERNATE_CODON:
                mCCV[var.mCodonIdx].mCounts[var.mValue] += 1;
                break;
            case VariationType::DELETION:
            case VariationType::INSERTION:
                {
                    bool isFrameShift = (var.mValue%3) != 0;
                    mCCV[var.mCodonIdx].mIndelCounts[isFrameShift] += 1;
                    if ( isFrameShift )
                        readHasFrameShiftIndel = true;
                    else
                        readHasNonFrameShiftIndel = true;
                }
                break;
            }
        }
        if ( readHasFrameShiftIndel )
            mNFSIndelReads += 1;
        if ( readHasNonFrameShiftIndel )
            mNFPIndelReads += 1;
    }

    void dumpNonAligner( char const* fileName, size_t readId,
                            Sequence const& seq, Quals const& quals )
    { if ( !mpJunkOS ) return;
      std::ostream& os = *mpJunkOS;
      os << '@' << fileName << ':' << readId << '\n';
      for ( Call call : seq ) os << charFor(call);
      os << "\n+\n";
      for ( Qual qual : quals ) os << (char)(qual+FASTQ_QUAL_OFFSET);
      os << '\n'; }

    Reference const& mRef;
    double mMinPerBaseSWScore;
    Qual mMinQ;
    unsigned mMaxSeqErrs;
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
    static Offset const SW_REF_PAD = 5;
};
Offset const ORFCaller::SW_REF_PAD;

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

void ORFCaller::writeIndelCounts( std::string const& filename )
{
    std::ofstream ofs(filename);
    ofs << "NFS\tFS\tTotCvg\n";

    for ( CodonCounts const& cCounts : mCCV )
    {
        ofs << cCounts.mIndelCounts[0] << '\t'
            << cCounts.mIndelCounts[1] << '\t'
            << cCounts.mCoverage << '\n';
    }
    ofs.close();
    if ( !ofs )
        FatalErr("Couldn't write indel counts to " << filename);
}

void ORFCaller::report( std::ostream& os )
{
    if ( !mNAlignedReads )
        FatalErr("Not even a single read aligned.");

    auto iORFWTCodon = mRef.getORFWTCodons().begin();
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
        size_t aaId = mAANamer.indexFor(*iORFWTCodon);
        ++iORFWTCodon;
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
        double dTot = cCounts.mCoverage/100.;
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

ReadReport ORFCaller::findIndels( Alignment const& aln,
                                  Sequence const& seq,
                                  Quals const& quals )
{
    Sequence const& ref = mRef.getRef();
    Offset refStart = aln.getRefStart()-std::min(aln.getRefStart(),SW_REF_PAD);
    Sequence::const_iterator refItr(ref.begin()+refStart);
    Sequence::const_iterator refEnd(ref.begin()+
                std::min(aln.getRefEnd()+SW_REF_PAD,ref.size()));
    Sequence::const_iterator seqItr, seqEnd;
    Quals::const_iterator qualItr;

    if ( !aln.isReadRC() )
    {
        seqItr = seq.begin() + aln.getReadStart();
        seqEnd = seq.begin() + aln.getReadEnd();
        qualItr = seq.begin() + aln.getReadStart();
    }
    else
    {
        mSeqTmp.clear();
        mSeqTmp.reserve(seq.size());
        for ( auto itr = seq.rbegin(), end = seq.rend(); itr != end; ++itr )
            mSeqTmp.push_back(complement(*itr));
        seqItr = mSeqTmp.begin() + aln.getReadStart();
        seqEnd = mSeqTmp.begin() + aln.getReadEnd();
        mQualsTmp = quals;
        std::reverse(mQualsTmp.begin(),mQualsTmp.end());
        qualItr = mQualsTmp.begin() + aln.getReadStart();
    }

    auto const& blocks = mSW.align(seqItr,seqEnd,refItr,refEnd);
    auto blkItr = blocks.begin();
    if ( blkItr->mDir == SmithWatermanizer::Trace::UP )
    {
        refStart += blkItr->mLen;
        refItr += blkItr->mLen;
        ++blkItr;
    }

    AlignmentScorer alignmentScorer(mRef,refStart,mMinQ);
    std::vector<Variation> indels;
    for ( auto blkEnd(blocks.end()); blkItr != blkEnd; ++blkItr )
    {
        unsigned nnn = blkItr->mLen;
        switch ( blkItr->mDir )
        {
        case SmithWatermanizer::Trace::UP:
            if ( alignmentScorer.inORF() &&
                    qualItr[-1] >= mMinQ && qualItr[0] >= mMinQ )
                indels.push_back(Variation{VariationType::DELETION,
                                            alignmentScorer.getCodonIdx(),
                                            nnn});
            while ( nnn-- )
            {
                alignmentScorer(*refItr,NCALL,0);
                ++refItr;
            }
            break;
        case SmithWatermanizer::Trace::LEFT:
            {
                bool loQ = false;
                while ( nnn-- )
                {
                    if ( *qualItr < mMinQ )
                        loQ = true;
                    ++seqItr; ++qualItr;
                }
                if ( !loQ && alignmentScorer.inORF() )
                    indels.push_back(Variation{VariationType::INSERTION,
                                                alignmentScorer.getCodonIdx(),
                                                blkItr->mLen});
            }
            break;
        case SmithWatermanizer::Trace::DIAG:
            while ( nnn-- )
            {
                alignmentScorer(*refItr,*seqItr,*qualItr);
                ++refItr; ++seqItr; ++qualItr;
            }
            break;
        case SmithWatermanizer::Trace::DONE:
            FatalErr("found DONE block in traceback");
        }
    }
    ReadReport report = alignmentScorer.getReport();
    report.mVariations.insert(report.mVariations.end(),
                                    indels.begin(),indels.end());
    return report;
}

void ORFCaller::processFastqs( char** beg, char** end )
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
                processRead(sq,*qItr,fileName,readId++);
                ++qItr;
            }
        }
        std::cerr << readId << " reads in file " << fileId << ".\n";
        fileId += 1;
    }
}

void ORFCaller::processPairedFastqs( char** beg, char** end )
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
                processPairedRead(*sItr1,*qItr1,*sItr2,*qItr2,
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
    unsigned maxSeqErrs = 2;
    double minPerBaseSWScore = 2.;
    std::string codonTranslation =
            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF";

    int opt;
    while ( (opt = getopt(argc, argv, "dE:K:M:O:pQ:S:T:")) != -1 )
    {
        switch ( opt )
        {
        case 'd':
            dumpNonAligners = true;
            break;
        case 'E':
        { std::istringstream iss(optarg);
          if ( !(iss >> maxSeqErrs) )
              FatalErr("Can't interpret E='" << optarg
                      << "'. Expecting an integer maximum for sequencing errors per read."); }
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
            if ( optopt == 'E' )
                FatalErr("Need value for E (max seq errs) option.");
            else if ( optopt == 'K' )
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

    if ( optind+2 > argc )
        FatalErr("Usage: ORFCall [-d] [-E MAX_SEQ_ERRS] [-K KMER_SIZE] "
                      " [-M MAX_MISMATCH_FRACTION]"
                      " [-O OUTPUT_FILENAME] [-p] [-Q MIN_QUAL] [-T CODON_TABLE] "
                      "ref.fasta codons_used.dat reads1.fastq [reads2.fastq...]\n"
                      "The -p flag signals paired mode.\n");

    Reference ref(argv[optind],argv[optind+1],K);
    ref.report();
    optind += 2;

    std::ofstream* pJunkOS = nullptr;
    if ( dumpNonAligners )
        pJunkOS = new std::ofstream(outputFilename+".unaligned.fasta");
    AminoAcidNamer aaNamer(codonTranslation);
    ORFCaller caller(ref,minPerBaseSWScore,minQ,maxSeqErrs,aaNamer,pJunkOS);
    char** beg = argv+optind;
    char** end = argv+argc;
    if ( !pairedMode )
        caller.processFastqs(beg,end);
    else
        caller.processPairedFastqs(beg,end);
    std::cerr << "Aligned reads: " << 100.*caller.getFractionAligned() << "%\n";

    if ( pJunkOS )
        delete pJunkOS;

    caller.writeCodonCounts(outputFilename+".cdn");
    caller.writeAminoAcidCounts(outputFilename+".aa");
    caller.writeIndelCounts(outputFilename+".indels");
    caller.report(std::cout);
}
