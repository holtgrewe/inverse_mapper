// ==========================================================================
//                             Inverse Mapper
// ==========================================================================
// Copyright (c) 2006-2012, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Smarter match management? Duplicate hits and unorderedly yield matches can pose problems. We could collect all matches of a fragment on a chromosome before any disabling.

// TODO(holtgrew): Allow -H to be a list of files.

#include <fstream>
#include <iostream>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/index/find_pigeonhole.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The Paths to the host genome.
    seqan::StringSet<seqan::CharString> hostGenomePaths;

    // The Path to the target genome.
    seqan::CharString targetGenomePath;

    // The word size.
    unsigned wordSize;

    // The minimal distance a read must have to not be filtered out.
    unsigned minDistance;
    // The maximal number of matches a read may have with its best distance before it is filtered out.
    unsigned maxBestHits;
    // The distance to use for the filtration.
    unsigned filtrationDistance;
    // The fraction of a fragment to consider as border, errors here do not count into XC tag.
    double borderFrac;
    // Whether or not to allow overlapping matches.
    bool allowOverlappingMatches;

    // The error rate, computed from word size and min distance.
    double errorRate;

    // Minimal length of repeats to mask.
    int repeatLength;

    // Path to output file.
    seqan::CharString outFile;

    AppOptions() : verbosity(1), wordSize(0), minDistance(0), maxBestHits(0), filtrationDistance(0), borderFrac(0.25), allowOverlappingMatches(false), errorRate(0), repeatLength(1000)
    {}

    unsigned computeQ() const
    {
        return wordSize / (minDistance + 1);
    }
};

// --------------------------------------------------------------------------
// Class ReadStats
// --------------------------------------------------------------------------

struct ReadStats
{
    // Sample position.
    int targetRefId;
    int targetPos;
    
    // Distance of best match, -1 means no match.
    int bestFoundDistance;

    // Number of matches with best distance.
    int numBestMatches;

    // Whether or not this read has been disabled.
    bool enabled;

    // Position of previous match.
    bool prevMatchRC;
    int prevMatchRefId;
    int prevMatchBeginPos;

    // Returns whether to append (1), to clear and append (2), or to do nothing (0).
    int update(int distance, int refId, bool rc, int beginPos, AppOptions const & options)
    {
        int res = 0;

        if (distance == bestFoundDistance)
        {
            if (!options.allowOverlappingMatches)
                if ((prevMatchRefId == refId) && (rc == prevMatchRC) && (beginPos - prevMatchBeginPos < (int)options.wordSize))
                    return 0;  // Skip, too close.
            res = 1;

            numBestMatches += 1;
            // Disable if too many hits and there is no space for improvement in distance.
            if (distance == (int)options.minDistance && numBestMatches > (int)options.maxBestHits)
                enabled = false;

            prevMatchRC = rc;
            prevMatchRefId = refId;
            prevMatchBeginPos = beginPos;
        }
        else if (bestFoundDistance == -1 || distance < bestFoundDistance)
        {
            res = (bestFoundDistance != -1) ? 2 : 1;
            bestFoundDistance = distance;
            numBestMatches = 1;

            if (distance < (int)options.minDistance)
                enabled = false;

            prevMatchRC = rc;
            prevMatchRefId = refId;
            prevMatchBeginPos = beginPos;
        }

        return res;
    }

    ReadStats() : targetRefId(-1), targetPos(-1), bestFoundDistance(-1), numBestMatches(0), enabled(true), prevMatchRC(false), prevMatchRefId(-1), prevMatchBeginPos(-1)
    {}

    ReadStats(int targetRefId, int targetPos) : targetRefId(targetRefId), targetPos(targetPos), bestFoundDistance(-1), numBestMatches(0), enabled(true), prevMatchRC(false), prevMatchRefId(-1), prevMatchBeginPos(-1)
    {}
};

// --------------------------------------------------------------------------
// Class MatchInfo
// --------------------------------------------------------------------------

// Store information about a match.

struct MatchInfo
{
    int readId;
    int distance;
    int refId;
    int beginPos;
    int endPos;
    seqan::String<seqan::CigarElement<> > cigarString;

    MatchInfo() : readId(-1), distance(-1), refId(-1), beginPos(-1), endPos(-1)
    {}

    MatchInfo(int _readId, int _distance, int _refId, int _beginPos, int _endPos) :
            readId(_readId), distance(_distance), refId(_refId), beginPos(_beginPos), endPos(_endPos)
    {}

    inline bool operator<(MatchInfo const & rhs) const
    {
        return (this->readId < rhs.readId) ||
                ((this->readId == rhs.readId) && (this->refId < rhs.refId)) ||
                ((this->readId == rhs.readId) && (this->refId == rhs.refId) && (this->beginPos < rhs.beginPos));
                
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function isFileNewerThan()
// ----------------------------------------------------------------------------

// Returns true iff both paths exists and the file at leftPath is newer than the file at rightPath.

#if PLATFORM_WINDOWS
#define STAT_FUNC _stat
#define STAT_STRUCT _stat
#else   // #if PLATFORM_WINDOWS
#define STAT_FUNC stat
#define STAT_STRUCT stat
#endif  // #if PLATFORM_WINDOWS

bool isFileNewerThan(char const * leftPath, char const * rightPath)
{
    struct STAT_STRUCT leftBuf;
    struct STAT_STRUCT rightBuf;

    if (STAT_FUNC(leftPath, &leftBuf) != 0)
        return false;
    if (STAT_FUNC(rightPath, &rightBuf) != 0)
        return false;

    return leftBuf.st_mtime > rightBuf.st_mtime;
}

#undef STAT_FUNC
#undef STAT_STRUCT

// ----------------------------------------------------------------------------
// Function getCigarString2()
// ----------------------------------------------------------------------------

namespace seqan {

// Uses = for equality and X for mismatch.

template <
    typename TCigar,
    typename TGaps1,
    typename TGaps2>
inline void
getCigarString2(
    TCigar &cigar,
    TGaps1 &gaps1,
    TGaps2 &gaps2,
    int splicedGapThresh = 20)
{
    typename Iterator<TGaps1>::Type it1 = iter(gaps1, 0);
    typename Iterator<TGaps2>::Type it2 = iter(gaps2, 0);
    clear(cigar);
    char op, lastOp = ' ';
    unsigned numOps = 0;

//  std::cout << gaps1 << std::endl;
//  std::cout << gaps2 << std::endl;
    for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
    {
        if (isGap(it1))
        {
            if (isGap(it2))
                op = 'P';
            else if (isClipped(it2))
                op = '?';
            else
                op = 'I';
        } 
        else if (isClipped(it1))
        {
            op = '?';
        }
        else 
        {
            if (isGap(it2))
                op = 'D';
            else if (isClipped(it2))
                op = 'S';
            else
                op = (*it1 == *it2) ? '=' : 'X';
        }
        
        // append CIGAR operation
        if (lastOp != op)
        {
            if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
                lastOp = 'N';
            if (numOps > 0)
                append(cigar, CigarElement<>(lastOp, numOps));
            numOps = 0;
            lastOp = op;
        }
        ++numOps;
    }
//  if (atEnd(it1) != atEnd(it2))
//        std::cerr << "Invalid pairwise alignment:" << std::endl << gaps1 << std::endl << gaps2 << std::endl;
    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
    if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
        lastOp = 'N';
    if (numOps > 0)
        append(cigar, CigarElement<>(lastOp, numOps));
}

}  // namespace seqan

// ----------------------------------------------------------------------------
// Function trimSeqHeaderToId()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Such a method is already in the SeqAn library, but in extras. Remove here when it is in core.
void trimSeqHeaderToId(seqan::CharString & header)
{
    unsigned i = 0;
    for (; i < length(header); ++i)
        if (isspace(header[i]))
            break;
    resize(header, i);
}

// --------------------------------------------------------------------------
// Function countCenterMatchInfo()
// --------------------------------------------------------------------------

// Count number of errors (non-'=') in the center of the read (within border of the begin/end of the string).

int countCenterMatchInfo(seqan::String<seqan::CigarElement<> > const & cigarString, double border = 0.25)
{
    unsigned len = 0;
    for (unsigned i = 0; i < length(cigarString); ++i)
    {
        if (cigarString[i].operation== 'M' || cigarString[i].operation== 'I' || cigarString[i].operation== '=' ||
            cigarString[i].operation== 'X')
            len += cigarString[i].count;
    }
    unsigned start = static_cast<unsigned>(ceil(border * len));
    unsigned stop = static_cast<unsigned>(floor((1.0 - border) * len));

    int res = 0;
    unsigned pos = 0;
    for (unsigned i = 0; i < length(cigarString); ++i)
    {
        if (pos >= start && pos < stop && cigarString[i].operation != '=')
            res += cigarString[i].count;
        if (cigarString[i].operation== 'M' || cigarString[i].operation== 'I' || cigarString[i].operation== '=' ||
            cigarString[i].operation== 'X')
            pos += cigarString[i].count;
    }

    return res;
}

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("inverse_mapper");
    // Set short description, version, and date.
    setShortDescription(parser, "Complement of read mapping.");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\n\fIOPTONS\\fP] \\fB-H\\fP \\fIHOST.fa\\fP \\fB-t\\fP \\fITARGET.fa\\fP");
    addDescription(parser, "Search for subsequences from \\fITARGET.fa\\fP of a configured size with minimal configured alignment distance towards \\fIHOST.fa\\fP.");
    addDescription(parser, "Any character 'U' (not case sensitive) in a host sequence is converted into a 'T'.  Any other character than 'C', 'G', 'A', 'T', 'U' (not case sensitive) is converted to 'N'");

    // Add Verbosity Options.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    addOption(parser, seqan::ArgParseOption("uv", "ultra-verbose", "Enable ultra verbose output."));
    hideOption(parser, "ultra-verbose");

    // Input / Output Options.
    addSection(parser, "Input / Output Options");
    addOption(parser, seqan::ArgParseOption("H", "host", "Host sequence.", seqan::ArgParseArgument::INPUTFILE, "HOST", true));
    setRequired(parser, "host");
    addOption(parser, seqan::ArgParseOption("t", "target", "Target sequence.", seqan::ArgParseArgument::INPUTFILE, "TARGET"));
    setRequired(parser, "target");
    addOption(parser, seqan::ArgParseOption("o", "out-file", "Path to output file.", seqan::ArgParseArgument::OUTPUTFILE, "FILE"));

    // Algorithm Options.
    addSection(parser, "Algorithm Options");
    addOption(parser, seqan::ArgParseOption("w", "word-size", "Word size.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "word-size", "18");
    addOption(parser, seqan::ArgParseOption("d", "min-distance", "Smallest edit distance a sequence of \\fITARGET\\fP must have to \\fIHOST\\fP.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "min-distance", "2");
    addOption(parser, seqan::ArgParseOption("f", "filtration-distance", "Distance to use for the string filter search.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "filtration-distance", "2");
    addOption(parser, seqan::ArgParseOption("m", "max-best-hits", "Maximal number of best hits before disabling.", seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "max-best-hits", "10");
    addOption(parser, seqan::ArgParseOption("", "border-frac", "Fraction of fragments to consider as border.  The number of border bases is rounded up in the algorithm.  Not counted into XC tag in output.", seqan::ArgParseArgument::DOUBLE, "REAL"));
    setDefaultValue(parser, "border-frac", "0.25");
    addOption(parser, seqan::ArgParseOption("", "overlapping-matches", "If this flag is set, overlapping matches are allowed."));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    if (isSet(parser, "ultra-verbose"))
        options.verbosity = 42;

    resize(options.hostGenomePaths, seqan::length(getOptionValues(parser, "host")));
    for (unsigned i = 0; i < length(options.hostGenomePaths); ++i)
        options.hostGenomePaths[i] = getOptionValues(parser, "host")[i];

    getOptionValue(options.targetGenomePath, parser, "target");
    getOptionValue(options.wordSize, parser, "word-size");
    getOptionValue(options.minDistance, parser, "min-distance");
    getOptionValue(options.filtrationDistance, parser, "filtration-distance");
    getOptionValue(options.maxBestHits, parser, "max-best-hits");
    options.errorRate = ((1.0 * options.minDistance) / options.wordSize) + 0.0001;
    getOptionValue(options.outFile, parser, "out-file");
    getOptionValue(options.borderFrac, parser, "border-frac");
    options.allowOverlappingMatches = isSet(parser, "overlapping-matches");

    if (options.computeQ() < 3)
    {
        std::cerr << "ERROR: q-gram size too slow.  Increase word size or number or "
                  << "errors. [q = floor(w / (min_distance + 1))].\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function containsN()
// --------------------------------------------------------------------------

// Return if sequence contains N.

template <typename TSequence>
bool containsN(TSequence const & seq)
{
    typedef typename seqan::Value<TSequence const>::Type TAlphabet;

    typedef typename seqan::Iterator<TSequence const, seqan::Rooted>::Type TIter;
    for (TIter it = begin(seq, seqan::Rooted()); !atEnd(it); goNext(it))
        if (*it == seqan::unknownValue<TAlphabet>())
            return true;

    return false;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    double programStartTime = sysTime();
    
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    if (options.verbosity > 0)
        std::cout << "INVERSE MAPPER\n"
                  << "==============\n\n";

    // Print the command line arguments back to the user.
    if (options.verbosity > 0)
    {
        std::cout << "__OPTIONS____________________________________________________________________\n"
                  << '\n'
                  << "VERBOSITY      \t" << options.verbosity << '\n'
                  << "HOST(S)        \t" << options.hostGenomePaths[0] << '\n';
        for (unsigned i = 1; i < length(options.hostGenomePaths); ++i)
            std::cout << "               \t" << options.hostGenomePaths[i] << '\n';
        std::cout << "TARGET         \t" << options.targetGenomePath << '\n'
                  << "WORD SIZE      \t" << options.wordSize << '\n'
                  << "FILTER DISTANCE\t" << options.filtrationDistance << '\n'
                  << "MIN DISTANCE   \t" << options.minDistance << '\n'
                  << "BORDER FRAC    \t" << options.borderFrac << '\n'
                  << "MAX BEST HITS  \t" << options.maxBestHits << '\n'
                  << "Q-GRAM SIZE    \t" << options.computeQ() << '\n'
                  << "ERROR RATE     \t[%] " << 100.0 * options.errorRate << "\n\n";
    }

    // -----------------------------------------------------------------------
    // Load / Build Host FAI Index.
    // -----------------------------------------------------------------------

    if (options.verbosity > 0)
        std::cout << "__START-UP I/O ______________________________________________________________\n\n";

    // Try to load index and create on the fly if necessary.
    seqan::String<seqan::FaiIndex> faiIndices;
    resize(faiIndices, length(options.hostGenomePaths));
    for (unsigned i = 0; i < length(options.hostGenomePaths); ++i)
    {
        if (options.verbosity > 0)
        {
            std::cout << "Genome          \t" << options.hostGenomePaths[i] << '\n'
                      << "Loading Index...\t" << std::flush;
        }
        bool rebuild = false;
        seqan::CharString faiPath = options.hostGenomePaths[i];
        append(faiPath, ".fai");
        if (isFileNewerThan(toCString(options.hostGenomePaths[i]), toCString(faiPath)))
        {
            rebuild = true;
            std::cout << "OUTDATED\n";
        }
        else if (seqan::read(faiIndices[i], toCString(options.hostGenomePaths[i])) != 0)
        {
            rebuild = true;
            std::cout << "NOT FOUND\n";
        }
        if (rebuild)
        {
            if (options.verbosity > 0)
                std::cout << "Building Index...\t" << std::flush;
            if (seqan::build(faiIndices[i], toCString(options.hostGenomePaths[i])) != 0)
            {
                if (options.verbosity > 0)
                    std::cerr << "\nERROR: Index could not be loaded or built.\n";
                return 1;
            }
            if (options.verbosity > 0)
                std::cout << "DONE\n"
                          << "Writing Index...\t";
            if (seqan::write(faiIndices[i]) != 0)
            {
                if (options.verbosity > 0)
                    std::cerr << "\nERROR: Index could not be written after building.\n";
                return 1;
            }
            if (options.verbosity > 0)
                std::cout << "DONE\n";
        }
        else
        {
            if (options.verbosity > 0)
                std::cout << "DONE\n";
        }
    }

    // -----------------------------------------------------------------------
    // Load Target Genome, Build Fragments.
    // -----------------------------------------------------------------------

    if (options.verbosity > 0)
        std::cerr << "\nBuilding fragments...\t";

    seqan::SequenceStream seqStream(toCString(options.targetGenomePath));
    seqan::StringSet<seqan::CharString> targetGenomeIds;
    seqan::StringSet<seqan::Dna5String> targetGenomeSeqs;
    // Genome Fragments and the sources of the fragments as (seq id, pos) pairs.
    seqan::StringSet<seqan::Pair<unsigned> > targetGenomeSources;
    seqan::StringSet<seqan::Dna5String> targetGenomeFragments;

    // Load the target genome fragment-wise.
    seqan::CharString id;
    seqan::Dna5String seq;
    seqan::String<ReadStats> readStats;  // TODO(holtgrew): Rename, because read == fragment.
    for (unsigned seqId = 0; !atEnd(seqStream); ++seqId)
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "\nERROR: Problem reading target sequence!\n";
            return 1;
        }

        trimSeqHeaderToId(id);
        appendValue(targetGenomeIds, id);
        appendValue(targetGenomeSeqs, seq);

        for (unsigned i = 0; i + options.wordSize - 1 < length(seq); ++i)
        {
            seqan::Infix<seqan::Dna5String>::Type targetInfix = infix(seq, i, i + options.wordSize);
            if (containsN(targetInfix))
                continue;  // Skip N's.
            appendValue(targetGenomeSources, seqan::Pair<unsigned>(seqId, i));
            appendValue(targetGenomeFragments, targetInfix);
            appendValue(readStats, ReadStats(seqId, i));
            back(readStats).enabled = true;
        }
    }

    if (options.verbosity > 0)
        std::cerr << "DONE\n\n"
                  << "FRAGMENTS\t" << length(targetGenomeFragments) << "\n\n"
                  << "__PERFORMING SEARCH__________________________________________________________\n\n";

    // We collect the best matches of each fragment.
    seqan::StringSet<seqan::String<MatchInfo> > matches;
    resize(matches, length(readStats));

    seqan::String<seqan::Pair<unsigned, unsigned> > refIdToFileIndexPair;  // map global ref id to <file, index in file>

    for (unsigned globalRefId = 0, fileId = 0; fileId < length(faiIndices); ++fileId)
    {
        // Run pigeonhole+banded Myers search for each contig.
        for (unsigned refId = 0; refId < numSeqs(faiIndices[fileId]); ++refId)
        {
            appendValue(refIdToFileIndexPair, seqan::Pair<unsigned>(fileId, refId));
            
            seqan::Dna5String contigSeq;
            reserve(contigSeq, sequenceLength(faiIndices[fileId], refId));
            // Load contig, converting U's to T's.  Do so in a block to free the buffer after assignment to contigSeq.
            {
                seqan::CharString buffer;
                if (readSequence(buffer, faiIndices[fileId], refId) != 0)
                {
                    std::cerr << "ERROR: Could not load sequence " << fileId << ":" << refId << "\n";
                    return 1;
                }
                for (unsigned i = 0; i < length(buffer); ++i)
                    if (buffer[i] == 'u' || buffer[i] == 'U')
                        buffer[i] = 'T';
                contigSeq = buffer;
            }

            typedef seqan::Finder<seqan::Dna5String, seqan::Pigeonhole<> > TFilterFinder;

            typedef seqan::StringSet<seqan::Dna5String, seqan::Dependent<> > TFragmentSet;
            typedef TFragmentSet TReadSet;
            typedef seqan::Value<TReadSet>::Type /*const*/ TReadSeq;
            typedef seqan::Value<TReadSet>::Type /*const*/ TReadPrefix;
            typedef seqan::ModifiedString<TReadPrefix, seqan::ModReverse> TRevReadPrefix;

            typedef seqan::Dna5String TContigSeq;
            typedef seqan::Infix<TContigSeq>::Type TGenomeInfix;
            typedef seqan::Position<TGenomeInfix>::Type			TPosition;
            typedef seqan::ModifiedString<TGenomeInfix, seqan::ModReverse> TGenomeInfixRev;
            typedef seqan::Finder<TGenomeInfix> TMyersFinder;

            typedef seqan::Finder<TGenomeInfix> TMyersFinder;
            typedef seqan::Finder<TGenomeInfixRev>							TMyersFinderRev;
            typedef seqan::PatternState_<TReadPrefix,	seqan::Myers<seqan::AlignTextBanded<seqan::FindInfix, seqan::NMatchesNone_, seqan::NMatchesNone_>, seqan::True, void> > TPatternState;
            typedef seqan::PatternState_<TRevReadPrefix, seqan::Myers<seqan::AlignTextBanded<seqan::FindPrefix, seqan::NMatchesNone_, seqan::NMatchesNone_>, seqan::True, void> > TRPatternState;

            TPatternState	patternState;
            TRPatternState  revPatternState;

            unsigned        rightClip = 0;
            unsigned        contigLength = length(contigSeq);

            // Perform filtration on forward and reverse strand.
            for (unsigned pass = 0; pass < 2; ++pass)  // Pass 1: FWD, Pass 2: REV.
            {
                unsigned q = options.computeQ();
                if (options.verbosity >= 3)
                    std::cerr << "Building q-gram index with q = " << q << "...\n";
                typedef seqan::Shape<seqan::Dna5, seqan::SimpleShape>    TShape;
                TShape shape(q);
                if (options.verbosity >= 3)
                {
                    seqan::CharString bitstr;
                    shapeToString(bitstr, shape);
                    std::cerr << "SHAPE\t" << bitstr << '\n';
                }

                TFragmentSet fragments;
                seqan::String<unsigned> fragIdMap;  // Position in fragments to position in targetGenomeFragments.
                for (unsigned fragId = 0; fragId < length(targetGenomeFragments); ++fragId)
                {
                    if (readStats[fragId].enabled)
                    {
                        appendValue(fragIdMap, fragId);
                        appendValue(fragments, targetGenomeFragments[fragId]);
                    }
                }
                if (options.verbosity >= 2)
                    std::cerr << "LENGTH FRAGMENTS\t" << length(fragments) << '\n';

                typedef seqan::IndexQGram<TShape, seqan::OpenAddressing>  TSpec;
                typedef seqan::Index<TFragmentSet, TSpec>                 TQGramIndex;
                TQGramIndex index(fragments, shape);

                typedef seqan::Pattern<TQGramIndex, seqan::Pigeonhole<> > TFilterPattern;
                TFilterPattern filterPattern(index);
                filterPattern.params.printDots = (options.verbosity > 0);
                _patternInit(filterPattern, options.errorRate);
                indexRequire(host(filterPattern), seqan::QGramSADir());

                bool forward = (pass == 0);
                if (pass == 1)
                    reverseComplement(contigSeq);
                TFilterFinder filterFinder(contigSeq, options.repeatLength, 1);

                double contigStartTime = sysTime();
                if (options.verbosity > 0)
                    std::cerr << '[' << faiIndices[fileId].indexEntryStore[refId].name << "/" << (pass ? "rev " : "fwd ") << length(fragments) << " " << (length(contigSeq) / 1000 / 1000) << "M] " << std::flush;

                // The Align object to store the semiglobal alignmetns in.
                seqan::Align<seqan::Dna5String> align;

                while (find(filterFinder, filterPattern, options.errorRate))
                {
                    unsigned readId = fragIdMap[position(filterPattern).i1];
                    if (options.verbosity >= 4)
                        std::cerr << "FOUND\n";
                    // Skip if disabled.
                    if (!readStats[readId].enabled)
                        continue;

                    // Compute clipping for banded Myers if SWIFT it is left/right of the text.
                    patternState.leftClip = (beginPosition(filterFinder) >= 0) ? 0 : -beginPosition(filterFinder);
                    rightClip = (endPosition(filterFinder) <= contigLength) ? 0 : endPosition(filterFinder) - contigLength;

                    // Verify SWIFT hit.
                    //
                    // We use a much simpler version of verification than in RazerS for now.

                    // First, search for the most promising end position.

                    unsigned ndlLength = length(targetGenomeFragments[readId]);
                    int minScore = -static_cast<int>(options.filtrationDistance);

                    TGenomeInfix inf = infix(filterFinder);  // TODO(holtgrew): Cannot pass to myersFinder directly because of const issues.
                    TGenomeInfix origInf(inf);
                    setEndPosition(inf, endPosition(inf) - 1);
                    ndlLength -= 1;
                    TReadPrefix readPrefix(targetGenomeFragments[readId], ndlLength);
                    TMyersFinder myersFinder(inf);
                    int bestScore = seqan::MinValue<int>::VALUE;
                    unsigned bestPos = 0;
                    if (options.verbosity >= 4)
                        std::cerr << "    CONTIG: " << origInf << '\n'
                                  << "    READ:   " << targetGenomeFragments[readId] << '\n'
                                  << "    PREFIX: " << readPrefix << '\n';

                    while (find(myersFinder, readPrefix, patternState, minScore))  // Iterate over all start positions.
                    {
                        TPosition const pos = position(hostIterator(myersFinder));
                        int score = getScore(patternState);

                        SEQAN_ASSERT_LT(pos + 1, length(origInf));
                        if (origInf[pos + 1] != back(targetGenomeFragments[readId]) || origInf[pos + 1] == seqan::unknownValue<seqan::Dna5>())
                            if (--score < minScore)
                                continue;

                        if (getScore(patternState) > bestScore)
                        {
                            bestScore = score;
                            bestPos = pos;
                        }
                        if (options.verbosity >= 4)
                            std::cerr << "MYERS CANDIDATE\t" << position(filterPattern).i1 << "\t" << fileId << ":" << refId << "\t" << (forward ? '+' : '-') << "\t" << contigLength << "\t" << (beginPosition(inf) + pos) << '\t' << score << std::endl;
                    }
                    if (bestScore == seqan::MinValue<int>::VALUE)
                        continue;  // No Myers hit, look for next SWIFT hit.

                    // Second, search for the leftmost start position of the hit.
                    __int64 infEndPos = endPosition(inf);
                    __int64 newInfEndPos = beginPosition(inf) + bestPos + 1;
                    revPatternState.leftClip = infEndPos - newInfEndPos + rightClip;
                    setEndPosition(inf, newInfEndPos);
                    if (endPosition(inf) > (unsigned)(ndlLength - bestScore))
                        setBeginPosition(inf, endPosition(inf) - ndlLength + bestScore);
                    else
                        setBeginPosition(inf, 0);

                    // Correct bestScore for reverse search below.
                    if (origInf[bestPos + 1] != back(targetGenomeFragments[readId]) || origInf[bestPos + 1] == seqan::unknownValue<seqan::Dna5>())
                        bestScore += 1;

                    TRevReadPrefix readRev(readPrefix);
                    TGenomeInfixRev infRev(inf);
                    TMyersFinderRev myersFinderRev(infRev);
                    __int64 beginPos = newInfEndPos;
                    if (options.verbosity >= 4)
                        std::cerr << "  R CONTIG:     " << infRev << std::endl
                                  << "  R READPRFX:   " << readRev << std::endl
                                  << "    best score: " << bestScore << std::endl;
                    while (find(myersFinderRev, readRev, revPatternState, bestScore))
                        beginPos = newInfEndPos - (position(myersFinderRev) + 1);
                    if (beginPos == newInfEndPos)
                        continue;  // No Banded Myers hit, skip.

                    // Correct bestScore for output.
                    if (origInf[bestPos + 1] != back(targetGenomeFragments[readId]) || origInf[bestPos + 1] == seqan::unknownValue<seqan::Dna5>())
                        bestScore -= 1;

                    // Update statistics for match.
                    int matchBeginPos = forward ? beginPos : contigLength - (newInfEndPos + 1);
                    int matchEndPos = forward ? newInfEndPos + 1 : contigLength - beginPos;
                    // int oldBestFoundDistance = readStats[readId].bestFoundDistance;
                    int what = readStats[readId].update(-bestScore, globalRefId, !forward, matchBeginPos, options);

                    // Clear matches if read was disabled or we have a new best distance.
                    if (what == 2)
                        clear(matches[readId]);
                    else if (what == 0)
                        continue;
                    // Perform and store alignment.
                    resize(rows(align), 2);
                    assignSource(row(align, 0), infix(contigSeq, beginPos, newInfEndPos + 1));
                    if (!forward)
                        std::swap(matchBeginPos, matchEndPos);
                    assignSource(row(align, 1), targetGenomeFragments[readId]);
                    seqan::Score<int> scoringScheme(0, -999, -1001, -1000);
                    seqan::AlignConfig<true, false, false, true> alignConfig;
                    if (options.verbosity >= 3)
                        std::cerr << "matchBeginPos == " << matchBeginPos << ", matchEndPos == " << matchEndPos << "\n";
                    int aliScore = globalAlignment(align, scoringScheme, alignConfig, seqan::Gotoh());
                    if (options.verbosity >= 3)
                    {
                        std::cerr << "aliScore == " << aliScore << "\n";
                        std::cerr << align;
                    }
                    // SEQAN_ASSERT_EQ(-aliScore/scoreMismatch(scoringScheme), bestScore);
                    // TODO(holtgrew): The assertion above should work!
                    // Record match for read.
                    appendValue(matches[readId], MatchInfo(readId, -bestScore, globalRefId, matchBeginPos, matchEndPos));
                    SEQAN_ASSERT_EQ((int)length(matches[readId]), readStats[readId].numBestMatches);
                    seqan::getCigarString2(back(matches[readId]).cigarString, row(align, 0), row(align, 1));
                    if (front(back(matches[readId]).cigarString).operation == 'D')
                        erase(back(matches[readId]).cigarString, 0);
                    if (back(back(matches[readId]).cigarString).operation == 'D')
                        eraseBack(back(matches[readId]).cigarString);

                    // Store the single-end match in the list of raw matches.
                    /*SingleEndMatch match;
                      match.readId = readId;
                      match.strand = forward ? '+' : '-';
                      match.contigId = contigId;
                      match.beginPos = forward ? beginPos : contigLength - (newInfEndPos + 1);
                      match.endPos = forward ? newInfEndPos + 1 : contigLength - beginPos;
                      match.score = bestScore;
                      appendValue(rawMatches, match);*/
                }
                double contigTime = sysTime() - contigStartTime;
                double megaBasesPerSecond = 1.0 * length(contigSeq) / 1000 / contigTime;
                if (options.verbosity > 0)
                {
                    char buffer[100];
                    sprintf(buffer, "%.2f", contigTime);
                    std::cerr << ' ' << contigTime << "s";
                    sprintf(buffer, "%.2f", megaBasesPerSecond);
                    std::cerr << ' ' << buffer << "kbp/s\n";
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Output
    // -----------------------------------------------------------------------

    if (options.verbosity > 0)
        std::cout << "\n__WRITING OUTPUT_____________________________________________________________\n\n";

    std::fstream outFileStream;
    std::ostream * outStream = &std::cout;
    if (!empty(options.outFile))
    {
        std::cout << "OUTPUT FILE\t" << options.outFile << "\n";
        outFileStream.open(toCString(options.outFile), std::ios_base::out | std::ios_base::binary);
        if (!outFileStream.good())
        {
            std::cerr << "ERROR: Could not open output file " << options.outFile << '\n';
            return false;
        }
        outStream = &outFileStream;
    }
    else
    {
        std::cout << "OUTPUT STDOUT\n";
    }


    // Build BAM Header.
    //
    seqan::BamHeader bamHeader;
    resize(bamHeader.records, 2 + length(refIdToFileIndexPair));
    // @HD VN:1.4
    seqan::BamHeaderRecord & headerHD = bamHeader.records[0];
    headerHD.type = seqan::BAM_HEADER_FIRST;
    resize(headerHD.tags, 1);
    headerHD.tags[0].i1 = "VN";
    headerHD.tags[0].i2 = "1.4";
    // @PG ID:initial_search PN:inverse_mapper VN:0.1 CL:inverse_mapper ...
    seqan::BamHeaderRecord & headerPG = bamHeader.records[1];
    headerPG.type = seqan::BAM_HEADER_PROGRAM;
    resize(headerPG.tags, 4);
    headerPG.tags[0].i1 = "ID";
    headerPG.tags[0].i2 = "initial_search";
    headerPG.tags[1].i1 = "PN";
    headerPG.tags[1].i2 = "inverse_mapper";
    headerPG.tags[2].i1 = "VN";
    headerPG.tags[2].i2 = "0.1";
    headerPG.tags[3].i1 = "CL";
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0)
            appendValue(headerPG.tags[3].i2, ' ');
        append(headerPG.tags[3].i2, argv[i]);
    }
    // @SQ headers and sequenceInfos.
    seqan::StringSet<seqan::CharString> refNameStore;
    resize(bamHeader.sequenceInfos, length(refIdToFileIndexPair));
    for (unsigned i = 0, idx = 0; i < length(faiIndices); ++i)
    {
        for (unsigned j = 0; j < numSeqs(faiIndices[i]); ++j, ++idx)
        {
            appendValue(refNameStore, sequenceName(faiIndices[i], j));

            bamHeader.sequenceInfos[idx].i1 = sequenceName(faiIndices[i], j);
            bamHeader.sequenceInfos[idx].i2 = sequenceLength(faiIndices[i], j);

            bamHeader.records[2 + idx].type = seqan::BAM_HEADER_REFERENCE;
            resize(bamHeader.records[2 + idx].tags, 2);
            bamHeader.records[2 + idx].tags[0].i1 = "SN";
            bamHeader.records[2 + idx].tags[0].i2 = sequenceName(faiIndices[i], j);
            bamHeader.records[2 + idx].tags[1].i1 = "LN";
            std::stringstream ss;
            ss << sequenceLength(faiIndices[i], j);
            bamHeader.records[2 + idx].tags[1].i2 = ss.str();
        }
    }
    seqan::NameStoreCache<seqan::StringSet<seqan::CharString> > refNameStoreCache(refNameStore);
    seqan::BamIOContext<seqan::StringSet<seqan::CharString> > bamIOContext(refNameStore, refNameStoreCache);

    seqan::write2(*outStream, bamHeader, bamIOContext, seqan::Sam());

    seqan::BamAlignmentRecord record;
    for (unsigned readId = 0; readId < length(readStats); ++readId)
    {
        clear(record);
        ReadStats const & stats = readStats[readId];

        // Fill parts of record that is the same for all matches of this fragment.
        //
        // QNAME
        std::stringstream ss;
        ss << targetGenomeIds[stats.targetRefId] << "_" << stats.targetPos << "_" << (stats.targetPos + options.wordSize);
        record.qName = ss.str();
        // SEQ
        record.seq = targetGenomeFragments[readId];

        // Write out not-found record.
        if (empty(matches[readId]))
        {
            record.flag = seqan::BAM_FLAG_UNMAPPED;
            write2(*outStream, record, bamIOContext, seqan::Sam());
        }

        // Sort matches for the read.
        std::sort(begin(matches[readId], seqan::Standard()), end(matches[readId], seqan::Standard()));

        // Write out matches.
        for (unsigned i = 0; i < length(matches[readId]); ++i)
        {
            MatchInfo const & matchInfo = matches[readId][i];
            // FLAG
            record.flag = (matchInfo.beginPos > matchInfo.endPos) ? seqan::BAM_FLAG_RC : 0;
            // CIGAR
            record.cigar = matchInfo.cigarString;
            // REF
            record.rId = matchInfo.refId;
            // POS
            record.pos = std::min(matchInfo.beginPos, matchInfo.endPos);
            // TAGS
            seqan::BamTagsDict tags(record.tags);
            // Distance to sequence.
            setTagValue(tags, "NM", matchInfo.distance);
            // Number of hits.
            setTagValue(tags, "NH", stats.numBestMatches);
            // Number of errors in the two center quarters of a read.
            setTagValue(tags, "XC", countCenterMatchInfo(record.cigar, options.borderFrac));
            // Whether or not the read was disabled.
            setTagValue(tags, "XD", (int)!stats.enabled);
            
            write2(*outStream, record, bamIOContext, seqan::Sam());
        }
        //if (readStats[i].bestFoundDistance != -1 && readStats[i].bestFoundDistance < (int)options.minDistance)
        //    continue;  // Skip.

        // (*outStream) << targetGenomeFragments[readId] << '\t' << targetGenomeIds[targetGenomeSources[readId].i1]
        //              << '\t' << targetGenomeSources[readId].i2 << '\t' << readStats[readId].bestFoundDistance
        //              << '\t' << readStats[readId].numBestMatches << '\n';
    }

    if (options.verbosity >= 1)
        std::cerr << "Total Time: " << (sysTime() - programStartTime) << " s\n";
    
    return 0;
}
