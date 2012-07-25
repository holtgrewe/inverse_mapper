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

#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>

#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/index/find_pigeonhole.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

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

    // The Path to the host genome.
    seqan::CharString hostGenomePath;

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

    // The error rate, computed from word size and min distance.
    double errorRate;

    // Minimal length of repeats to mask.
    int repeatLength;

    // Path to output file.
    seqan::CharString outFile;

    AppOptions() : verbosity(1), wordSize(0), minDistance(0), maxBestHits(0), filtrationDistance(0), errorRate(0), repeatLength(1000)
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

    void update(int distance, int refId, bool rc, int beginPos, AppOptions const & options)
    {
        if (distance == bestFoundDistance)
        {
            if ((prevMatchRefId == refId) && (rc == prevMatchRC) && (beginPos - prevMatchBeginPos < (int)options.wordSize))
                return;  // Skip, too near.

            numBestMatches += 1;
            // TODO(holtgrew): We should only purge these ones if there is no room for improvement, right?
            if (numBestMatches > (int)options.maxBestHits)
                enabled = false;

            prevMatchRC = rc;
            prevMatchRefId = refId;
            prevMatchBeginPos = beginPos;
        }
        else if (bestFoundDistance == -1 || distance < bestFoundDistance)
        {
            bestFoundDistance = distance;
            numBestMatches = 1;

            if (distance < (int)options.minDistance)
                enabled = false;

            prevMatchRC = rc;
            prevMatchRefId = refId;
            prevMatchBeginPos = beginPos;
        }
    }

    ReadStats() : bestFoundDistance(-1), numBestMatches(0), enabled(true), prevMatchRC(false), prevMatchRefId(-1), prevMatchBeginPos(-1)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

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
    addUsageLine(parser, "[\\fIOPTIONS\\fP] [\\fIOPTONS\\fP] \\fB-H\\fP \\fIHOST.fa\\fP \\fB-V\\fP \\fITARGET.fa\\fP");
    addDescription(parser, "Search for subsequences from \\fITARGET.fa\\fP of a configured size with minimal configured alignment distance towards \\fIHOST.fa\\fP.");

    // Add Verbosity Options.
    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    // Input / Output Options.
    addSection(parser, "Input / Output Options");
    addOption(parser, seqan::ArgParseOption("H", "host", "Host sequence.", seqan::ArgParseArgument::INPUTFILE, "HOST"));
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

    getOptionValue(options.hostGenomePath, parser, "host");
    getOptionValue(options.targetGenomePath, parser, "target");
    getOptionValue(options.wordSize, parser, "word-size");
    getOptionValue(options.minDistance, parser, "min-distance");
    getOptionValue(options.filtrationDistance, parser, "filtration-distance");
    getOptionValue(options.maxBestHits, parser, "max-best-hits");
    options.errorRate = ((1.0 * options.minDistance) / options.wordSize) + 0.0001;
    getOptionValue(options.outFile, parser, "out-file");

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
                  << "HOST           \t" << options.hostGenomePath << '\n'
                  << "TARGET         \t" << options.targetGenomePath << '\n'
                  << "WORD SIZE      \t" << options.wordSize << '\n'
                  << "FILTER DISTANCE\t" << options.filtrationDistance << '\n'
                  << "MIN DISTANCE   \t" << options.minDistance << '\n'
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
    seqan::FaiIndex faiIndex;
    if (options.verbosity > 0)
        std::cout << "Loading Index...\t" << std::flush;
    if (seqan::read(faiIndex, toCString(options.hostGenomePath)) != 0)
    {
        if (options.verbosity > 0)
            std::cout << "NOT FOUND\n"
                      << "Building Index...\t" << std::flush;
        if (seqan::build(faiIndex, toCString(options.hostGenomePath)) != 0)
        {
            if (options.verbosity > 0)
                std::cerr << "\nERROR: Index could not be loaded or built.\n";
            return 1;
        }
        if (options.verbosity > 0)
            std::cout << "DONE\n"
                      << "Writing Index...\t";
        if (seqan::write(faiIndex) != 0)
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

    // -----------------------------------------------------------------------
    // Load Target Genome, Build Fragments.
    // -----------------------------------------------------------------------

    if (options.verbosity > 0)
        std::cerr << "Building fragments...\t";

    seqan::SequenceStream seqStream(options.targetGenomePath);
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
            appendValue(readStats, ReadStats());
        }
    }

    if (options.verbosity > 0)
        std::cerr << "DONE\n\n"
                  << "FRAGMENTS\t" << length(targetGenomeFragments) << "\n\n"
                  << "__PERFORMING SEARCH__________________________________________________________\n\n";

    // Run pigeonhole+banded Myers search for each contig.
    for (unsigned refId = 0; refId < numSeqs(faiIndex); ++refId)
    {
        // Load contig.
        seqan::Dna5String contigSeq;
        if (readSequence(contigSeq, faiIndex, refId) != 0)
        {
            std::cerr << "ERROR: Could not load sequence " << refId << "\n";
            return 1;
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
            // std::cerr << "Building q-gram index with q = " << q << "...\n";
            typedef seqan::Shape<seqan::Dna5, seqan::SimpleShape>    TShape;
            TShape shape(q);

            TFragmentSet fragments;
            for (unsigned fragId = 0; fragId < length(targetGenomeFragments); ++fragId)
                if (readStats[fragId].enabled)
                    assignValueById(fragments, targetGenomeFragments[fragId]);

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
                std::cerr << '[' << faiIndex.indexEntryStore[refId].name << "/" << (pass ? "rev " : "fwd ") << length(fragments) << " " << (length(contigSeq) / 1000 / 1000) << "M] " << std::flush;
            while (find(filterFinder, filterPattern, options.errorRate))
            {
                unsigned readId = idToPosition(targetGenomeFragments, positionToId(fragments, position(filterPattern).i1));
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

                unsigned ndlLength = sequenceLength(readId, targetGenomeFragments);
                int minScore = -static_cast<int>(options.filtrationDistance);

                TGenomeInfix inf = infix(filterFinder);  // TODO(holtgrew): Cannot pass to myersFinder directly because of const issues.
                TGenomeInfix origInf(inf);
                setEndPosition(inf, endPosition(inf) - 1);
                ndlLength -= 1;
                TReadPrefix readPrefix(targetGenomeFragments[readId], ndlLength);
                TMyersFinder myersFinder(inf);
                int bestScore = seqan::MinValue<int>::VALUE;
                unsigned bestPos = 0;
                if (options.verbosity >= 11)
                    std::cerr << "    CONTIG: " << origInf << std::endl
                              << "    READ:   " << targetGenomeFragments[readId] << std::endl
                              << "    PREFIX: " << readPrefix << std::endl;

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
                    if (options.verbosity >= 11)
                        std::cerr << "MYERS CANDIDATE\t" << position(filterPattern).i1 << "\t" << refId << "\t" << (forward ? '+' : '-') << "\t" << contigLength << "\t" << (beginPosition(inf) + pos) << std::endl;
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
                int pos = forward ? beginPos : contigLength - (newInfEndPos + 1);
                readStats[readId].update(-bestScore, refId, !forward, pos, options);

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

    (*outStream) << "fragment\ttarget ref\ttarget pos\tbest distance\tnum matches\n";
    for (unsigned i = 0; i < length(readStats); ++i)
    {
        //if (readStats[i].bestFoundDistance != -1 && readStats[i].bestFoundDistance < (int)options.minDistance)
        //    continue;  // Skip.

        (*outStream) << targetGenomeFragments[i] << '\t' << targetGenomeIds[targetGenomeSources[i].i1]
                     << '\t' << targetGenomeSources[i].i2 << '\t' << readStats[i].bestFoundDistance
                     << '\t' << readStats[i].numBestMatches << '\n';
    }

    return 0;
}
