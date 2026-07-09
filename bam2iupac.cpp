/* MIT License
 * Copyright (c) 2024 Kristian Ullrich
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software i
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <vector>
#include <getopt.h>
#include <cmath>
#include <tuple>
#include <fstream>
#include "sam.h"

#define VERSION "0.0.3"
#define EXENAME "bam2iupac"
#define GITHUB_URL "https://github.com/kullrich/bam2iupac"

#undef BAM_CIGAR_STR
#define BAM_CIGAR_STR "MIDNSHP=XB"

struct BamExtractionResult {
    std::vector<std::vector<int>> counts;
    std::string chromosome;
    hts_pos_t startPos;
    hts_pos_t endPos;
};

std::tuple<std::string, hts_pos_t, hts_pos_t> resolveRegion(
    bam_hdr_t* header,
    const std::string& chromosome,
    hts_pos_t startPos,
    hts_pos_t endPos) {

    int tid = -1;

    for (int i = 0; i < header->n_targets; ++i) {
        if (chromosome == header->target_name[i]) {
            tid = i;
            break;
        }
    }

    if (tid < 0) {
        throw std::runtime_error(
            "Chromosome not found in BAM header: " + chromosome);
    }

    hts_pos_t chromosomeLength = header->target_len[tid];

    // Special case: entire chromosome requested
    if (startPos == 0 && endPos == 0) {
        startPos = 1;
        endPos = chromosomeLength;
    } else {

        // Reject invalid start
        if (startPos < 1) {
            throw std::runtime_error(
                "Start position must be >= 1.");
        }

        // Clamp end to chromosome length
        if (endPos > chromosomeLength) {
            endPos = chromosomeLength;
        }

        // After clamping, the region may have become invalid
        if (startPos > endPos) {
            throw std::runtime_error(
                "Requested region starts beyond the end of chromosome.");
        }
    }

    return std::make_tuple(chromosome, startPos, endPos);
}

std::tuple<std::string, hts_pos_t, hts_pos_t> parseRegion(
    const std::string& regionStr) {

    size_t colonPos = regionStr.find(':');
    size_t dashPos = regionStr.find('-');
    size_t spacePos = regionStr.find(' ');
    std::string chromosome = "";
    hts_pos_t startPos = 0;
    hts_pos_t endPos = 0;

    if (colonPos != std::string::npos && dashPos != std::string::npos) {
        // Format is chr:start-end
        chromosome = regionStr.substr(static_cast<size_t>(0), colonPos);
        startPos = std::stoi(regionStr.substr(colonPos + static_cast<size_t>(1), dashPos - colonPos - static_cast<size_t>(1)));
        endPos = std::stoi(regionStr.substr(dashPos + static_cast<size_t>(1)));
    } else if (spacePos != std::string::npos) {
        // Format is chr start end
        std::istringstream iss(regionStr);
        iss >> chromosome >> startPos >> endPos;
    } else {
        throw std::runtime_error("Error: Invalid region format. Please use one of the following formats: 'chr:start-end' or 'chr start end'.");
    }

    if (!(startPos == 0 && endPos == 0)) {
        if (startPos < 1) {
            throw std::runtime_error("Error: Start position must be >= 1 (BAM coordinates are 1-based).");
        }
        if (startPos > endPos) {
            throw std::runtime_error("Error: Either 'startPos > endPos' or region not parsed correctly.");
        }
    }

    return std::make_tuple(chromosome, startPos, endPos);
}

//char intToIupac[15] = {
//    'A', // Adenine
//    'C', // Cytosine
//    'G', // Guanine
//    'T', // Thymine
//    'R', // puRine A+G
//    'Y', // pYrimidine C+T
//    'S', // Strong interaction (3 H bonds) G+C
//    'W', // Weak interaction (2 H bonds) A+T
//    'K', // Keto G+T
//    'M', // aMino A+C
//    'B', // not-A, B follows A
//    'D', // not-C, D follows C
//    'H', // not-G, H follows G in the alphabet
//    'V', // not-T (not-U), V follows U
//    'N'}; // aNy

static const char maskToIupac[16] = {
    'N', // 0000 (invalid/no base)
    'A', // 0001
    'C', // 0010
    'M', // 0011 A+C
    'G', // 0100
    'R', // 0101 A+G
    'S', // 0110 C+G
    'V', // 0111 A+C+G
    'T', // 1000
    'W', // 1001 A+T
    'Y', // 1010 C+T
    'H', // 1011 A+C+T
    'K', // 1100 G+T
    'D', // 1101 A+G+T
    'B', // 1110 C+G+T
    'N'  // 1111 any base
};

//char getIupac(const std::vector<int>& counts, double iRatio) {
//
//    int depth = 0;
//    int whichIUPAC = 0;
//    double bIUPACscore = 0.0;
//
//    for (int b = 0; b < 4; b++) {
//        depth += counts[static_cast<unsigned int>(b)];
//    }
//    if (depth <= 0) {
//        whichIUPAC = 14;
//    } else {
//        for (int b = 0; b < 4; b++) {
//            if (static_cast<double>(counts[static_cast<unsigned int>(b)]) / static_cast<double>(depth) > iRatio) {
//                bIUPACscore += pow(b + 1, 2);
//            }
//        }
//        //N;A;C;G;T;A+G;C+T;G+C;A+T;G+T;A+C;C+G+T;A+G+T;A+C+T;A+C+G;A+C+G+T
//        const double scoreThresholds[] = {0.0, 1.0, 4.0, 9.0, 16.0, 10.0, 20.0, 13.0, 17.0, 25.0, 5.0, 29.0, 26.0, 21.0, 14.0, 30.0};
//        const int iupacValues[] = {14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
//
//        for (size_t i = static_cast<size_t>(0); i < sizeof(scoreThresholds) / sizeof(scoreThresholds[0]); i++) {
//            if (bIUPACscore == scoreThresholds[i]) {
//                whichIUPAC = iupacValues[i];
//                break;
//            }
//        }
//    }
//    return intToIupac[whichIUPAC];
//}

char getIupac(
    const std::vector<int>& counts,
    double iRatio) {

    int depth = 0;

    // counts:
    // 0=A
    // 1=C
    // 2=G
    // 3=T
    // 4=N
    // 5=deletion

    for (int b = 0; b < 4; b++) {
        depth += counts[static_cast<size_t>(b)];
    }

    if (depth <= 0) {
        return 'N';
    }

    uint8_t mask = 0;
    double threshold = static_cast<double>(depth) * iRatio;
    if (counts[0] > threshold)
        mask |= 1;   // A
    if (counts[1] > threshold)
        mask |= 2;   // C
    if (counts[2] > threshold)
        mask |= 4;   // G
    if (counts[3] > threshold)
        mask |= 8;   // T

    return maskToIupac[mask];
}

std::string generateFasta(
    const std::vector<std::vector<std::vector<std::string>>>& allIupacStrings,
    const std::vector<std::string>& sequenceIds,
    const std::vector<std::string>& chromosomes,
    const std::vector<int>& startPositions,
    const std::vector<int>& endPositions) {

    std::stringstream fastaStream;

    for (size_t i = static_cast<size_t>(0); i < allIupacStrings[static_cast<size_t>(0)].size(); ++i) {

        const std::string &sequenceId = sequenceIds[i];
        std::string sequence;
        fastaStream << ">" << sequenceId << " ";

        for (size_t r = static_cast<size_t>(0); r < allIupacStrings.size(); ++r) {
            const std::string &regionInfo = chromosomes[r] + " " + std::to_string(static_cast<long long>(startPositions[r])) + " " + std::to_string(static_cast<long long>(endPositions[r]));
            fastaStream << regionInfo << " ";

            for (const std::string& nucleotide : allIupacStrings[r][i]) {
                sequence += nucleotide;
            }
        }
        fastaStream << std::endl;
        fastaStream << sequence << std::endl;
        sequence.clear();
    }

    return fastaStream.str();
}

void printIndCounts(
    const std::vector<std::vector<int>>& Counts) {

    for (size_t p = static_cast<size_t>(0); p < Counts.size(); ++p) {
        std::cout << "    Position " << p << ": ";
        for (size_t c = static_cast<size_t>(0); c < Counts[p].size(); ++c) {
            std::cout << Counts[p][c] << " ";
        }
        std::cout << "\n";
    }
}

void printAllCounts(
    const std::vector<std::string> bamFiles,
    const std::vector<std::string> chromosomes,
    const std::vector<int> startPositions,
    const std::vector<int> endPositions,
    const std::vector<std::vector<std::vector<std::vector<int>>>>& allCounts) {

    for (size_t r = static_cast<size_t>(0); r < allCounts.size(); ++r) {
        std::cout << "Region:" << "\n";
        std::cout << "\t" << chromosomes[r] << ":" << startPositions[r] << "-" << endPositions[r] << "\n";
        for (size_t b = static_cast<size_t>(0); b < allCounts[r].size(); ++b) {
            std::cout << "BAM file:" << "\n";
            std::cout << "\t" << bamFiles[b] << "\n";
            std::cout << "CHR\tPOS\tA\tC\tG\tT\tN\t-" << "\n";
            for (size_t p = static_cast<size_t>(0); p < allCounts[r][b].size(); ++p) {
                std::cout << chromosomes[r] << "\t" << static_cast<size_t>(startPositions[r]) + p << "\t";
                for (size_t c = static_cast<size_t>(0); c < allCounts[r][b][p].size(); ++c) {
                    std::cout << allCounts[r][b][p][c] << "\t";
                }
                std::cout << "\n";
            }
        }
    }
}

BamExtractionResult extractCountsFromBam(
    const std::string inputBamFile,
    std::string chromosome,
    hts_pos_t startPos,
    hts_pos_t endPos,
    const int minMapQuality,
    const int maxMapQuality,
    const int minBaseQuality,
    const int maxBaseQuality,
    const int minCoverage,
    const int maxCoverage,
    const int debug) {

    // Set nullptr
    samFile *bam = nullptr;
    bam_hdr_t *header = nullptr;
    hts_idx_t *idx = nullptr;
    hts_itr_t *iter = nullptr;

    // Initialize counts vector
    std::vector<std::vector<int>> counts; // Stores counts for each position
    // Initialize resolved region
    std::tuple<std::string, hts_pos_t, hts_pos_t> regionResolved;

    // Process BAM file
    try {
        // Open BAM file
        bam = sam_open(inputBamFile.c_str(), "r");
        if (bam == nullptr) {
            throw std::runtime_error("Failed to open BAM file.");
        }

        // Read header
        header = sam_hdr_read(bam);
        if (header == nullptr) {
            throw std::runtime_error("Failed to read BAM header.");
        }
        regionResolved = resolveRegion(header, chromosome, startPos, endPos);
        // overwrite coordinates with resolved values
        chromosome = std::get<0>(regionResolved);
        startPos = std::get<1>(regionResolved);
        endPos = std::get<2>(regionResolved);

        if (debug == 1) {
            std::cout << "Resolved region: "
                      << chromosome << ":"
                      << startPos << "-"
                      << endPos << std::endl;
        }
        size_t numCounts = static_cast<size_t>(endPos) - static_cast<size_t>(startPos) + static_cast<size_t>(1); // Get length from start to end
        // Initialize counts from start to end
        // A, C, G, T, N, -
        counts.resize(numCounts, std::vector<int>(static_cast<size_t>(6), 0));

        // Initialize BAM index
        idx = sam_index_load(bam, inputBamFile.c_str());
        if (idx == nullptr) {
            throw std::runtime_error("Failed to load BAM index.");
        }

        // Create iterator for the specified region
        iter = sam_itr_querys(
            idx, header, (chromosome + ":" + std::to_string(startPos) + "-" + std::to_string(endPos)).c_str());
        if (iter == nullptr) {
            throw std::runtime_error("Failed to create iterator for region.");
        }

        // Allocate memory for alignment
        bam1_t *alignment = bam_init1();

        // Process alignments within the specified region
        while (sam_itr_next(bam, iter, alignment) >= 0) {
            // Get details of the alignment
            const char *qname = bam_get_qname(alignment); // Get query name
            hts_pos_t pos = alignment->core.pos; // Position of alignment
            int32_t seqLen = alignment->core.l_qseq; // Length of the read sequence
            uint8_t mapQuality = alignment->core.qual; // Mapping quality of the alignment
            // uint16_t flag = alignment->core.flag; // Bitwise flag
            uint32_t n_cigar = alignment->core.n_cigar; // Number of CIGAR operations
            // int tid = alignment->core.tid; // ID of the reference sequence
            // int mtid = alignment->core.mtid; // ID of the mate reference sequence
            // int mpos = alignment->core.mpos; // Position of the mate
            // int isize = alignment->core.isize; // Insert size

            // Filter by mapping quality
            if (static_cast<int>(mapQuality) < minMapQuality || static_cast<int>(mapQuality) > maxMapQuality) {
                continue;
            }

            std::string algString; // Define the alignment sequence
            uint8_t *seq = bam_get_seq(alignment); // Get the sequence
            uint8_t *qual = bam_get_qual(alignment); // Get the base qualities
            uint32_t *cigar = bam_get_cigar(alignment); // Get the CIGAR operations

            if (cigar == nullptr) {
                std::cerr << "Failed to read CIGAR string." << std::endl;
            }

            /* bam_cigar_type returns a bit flag with:
             *   bit 1 set if the cigar operation consumes the query
             *   bit 2 set if the cigar operation consumes the reference
             *
             * For reference, the unobfuscated truth table for this function is:
             * BAM_CIGAR_TYPE  QUERY  REFERENCE STRING
             * --------------------------------
             * BAM_CMATCH      1      1      M
             * BAM_CINS        1      0      I
             * BAM_CDEL        0      1      D
             * BAM_CREF_SKIP   0      1      N
             * BAM_CSOFT_CLIP  1      0      S
             * BAM_CHARD_CLIP  0      0      H
             * BAM_CPAD        0      0      P
             * BAM_CEQUAL      1      1      =
             * BAM_CDIFF       1      1      X
             * BAM_CBACK       0      0      B
             * --------------------------------
             */

            hts_pos_t refPos = pos + 1;
            hts_pos_t seqPos = 0;

            if (debug == 1) {
                // Print the details of the alignment
                std::cout << "Query Name: " << qname << std::endl;
                std::cout << "Alignment Position: " << pos + 1 << "\t";
                std::cout << "Sequence Length: " << seqLen << "\t";
                std::cout << "Mapping Quality: " << mapQuality << "\t";
                std::cout << "Number of CIGAR Operations: " << n_cigar << std::endl;
                std::string cigarStr;
                for (uint32_t i = static_cast<uint32_t>(0); i < alignment->core.n_cigar; ++i) {
                    int opLen = bam_cigar_oplen(cigar[i]);
                    char op = bam_cigar_opchr(cigar[i]);
                    cigarStr += std::to_string(opLen) + op;
                }
                std::cout << "CIGAR: " << cigarStr << std::endl;
            }

            for (uint32_t i = static_cast<uint32_t>(0); i < alignment->core.n_cigar; ++i) {
                if (refPos > endPos) {
                    if (debug==1) {
                        std::cout << "Break since refPos (" << refPos << ") > endPos (" << endPos <<")." << std::endl;
                    }
                    break;
                }
                int opLen = bam_cigar_oplen(cigar[i]);
                int op = bam_cigar_op(cigar[i]);
                switch (op) {
                    case BAM_CMATCH:
                    case BAM_CEQUAL:
                    case BAM_CDIFF:
                        if (debug == 1) {
                            std::cout << "Query and Reference both changed" << std::endl;
                            std::cout << "opLen: " << opLen << std::endl;
                        }

                        for (int j = 0; j < opLen; ++j) {
                            char base = seq_nt16_str[bam_seqi(seq, seqPos)]; // Get base
                            int baseQuality = static_cast<int>(qual[seqPos]); // Get baseQuality
                            int baseIdx = refPos - startPos; // Get baseIdx
                            size_t bIdx = static_cast<size_t>(baseIdx);
                            algString += base;
                            if (refPos >= startPos && refPos <= endPos) {
                                if (baseQuality >= minBaseQuality && baseQuality <= maxBaseQuality) {
                                    switch (base) {
                                        case 'A': ++counts[bIdx][static_cast<size_t>(0)]; break;
                                        case 'C': ++counts[bIdx][static_cast<size_t>(1)]; break;
                                        case 'G': ++counts[bIdx][static_cast<size_t>(2)]; break;
                                        case 'T': ++counts[bIdx][static_cast<size_t>(3)]; break;
                                        case 'N': ++counts[bIdx][static_cast<size_t>(4)]; break;
                                        default: break;
                                    }
                                }
                            }
                            if (debug == 1) {
                                std::cout << "refPos: " << refPos << "\t";
                                std::cout << "seqPos: " << seqPos + 1 << "\t";
                                std::cout << "baseIdx: " << baseIdx << "\t";
                                std::cout << "base: " << base << "\t";
                                std::cout << "baseQ: " << baseQuality << std::endl;
                                //printIndCounts(counts);
                            }
                            refPos++;
                            seqPos++;
                            if (refPos > endPos) {
                                if (debug==1) {
                                    std::cout << "Break since refPos (" << refPos << ") > endPos (" << endPos <<")." << std::endl;
                                }
                                break;
                            }
                        }
                        break;
                    case BAM_CDEL:
                    case BAM_CREF_SKIP:
                        if (debug == 1) {
                            std::cout << "Only Reference changed" << std::endl;
                            std::cout << "opLen: " << opLen << std::endl;
                        }

                        for (int j = 0; j < opLen; ++j) {
                            int baseIdx = refPos - startPos; // Get baseIdx
                            size_t bIdx = static_cast<size_t>(baseIdx);
                            algString += "-";
                            if (refPos >= startPos && refPos <= endPos) {
                                // Increase Deletion count
                                ++counts[bIdx][static_cast<size_t>(5)];
                            }
                            if (debug == 1) {
                                std::cout << "refPos: " << refPos << "\t";
                                std::cout << "seqPos: " << seqPos + 1 << "\t";
                                std::cout << "baseIdx: " << baseIdx << "\t";
                                std::cout << "base: - " << "\t";
                                std::cout << "baseQ: *" << std::endl;
                                //printIndCounts(counts);
                            }
                            refPos++;
                            if (refPos > endPos) {
                                if (debug==1) {
                                    std::cout << "Break since refPos (" << refPos << ") > endPos (" << endPos <<")." << std::endl;
                                }
                                break;
                            }
                        }
                        break;
                    case BAM_CINS:
                    case BAM_CSOFT_CLIP:
                        if (debug == 1) {
                            std::cout << "Only Query changed" << std::endl;
                            std::cout << "opLen: " << opLen << std::endl;
                        }

                        for (int j = 0; j < opLen; ++j) {
                            algString += tolower(seq_nt16_str[bam_seqi(seq, seqPos )]);
                            seqPos++;
                        }
                        break;
                    case BAM_CHARD_CLIP:
                    case BAM_CPAD:
                    case BAM_CBACK:
                        if (debug == 1) {
                            std::cout << "Neither Query nor Reference changed" << std::endl;
                        }
                        break;
                    case 10:
                    case 11:
                    case 12:
                    case 13:
                    case 14:
                    case 15:
                        if (debug == 1) {
                            std::cout << "Reserved or user-defined operation: " << op << std::endl;
                        }
                        break;
                    default:
                        std::cerr << "Unsupported CIGAR operation: " << op << std::endl;
                        break;
                }
            }
                if (debug == 1) {
                    std::string seqString;
                    uint8_t *seq = bam_get_seq(alignment); // Get the sequence
                    for (int i = 0; i < seqLen; ++i) {
                        char base = seq_nt16_str[bam_seqi(seq, i)]; // Get base
                        seqString.push_back(base);
                    }
                    std::cout << "Sequence: " << std::endl;
                    std::cout << seqString << std::endl;
                    std::cout << "Aligned Sequence: " << std::endl;
                    std::cout << algString << std::endl << std::endl;
                }

        }

        // Cleanup
        bam_destroy1(alignment);
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(bam);
    } catch (const std::runtime_error& e) {
        //std::cerr << e.what() << std::endl;

        // Clean up resources in case of error
        if (iter != nullptr) {
            hts_itr_destroy(iter);
        }
        if (idx != nullptr) {
            hts_idx_destroy(idx);
        }
        if (header != nullptr) {
            bam_hdr_destroy(header);
        }
        if (bam != nullptr) {
            sam_close(bam);
        }

        // Optionally rethrow the exception or handle it
        throw;
    }

    // Apply minCoverage and maxCoverage filter
    for (size_t i = static_cast<size_t>(0); i < counts.size(); ++i) {
        int baseCoverage = counts[i][0] + counts[i][1] + counts[i][2] + counts[i][3] + counts[i][4];
        if (baseCoverage < minCoverage || baseCoverage > maxCoverage) {
            counts[i][0] = 0; // Set A count to 0
            counts[i][1] = 0; // Set C count to 0
            counts[i][2] = 0; // Set G count to 0
            counts[i][3] = 0; // Set T count to 0
            counts[i][4] = 0; // Set N count to 0
        }
    }

    return {
        counts,
        chromosome,
        startPos,
        endPos
    };
}

void show_help(
    const char* program_name,
    int retcode) {

    FILE* out = (retcode == EXIT_SUCCESS ? stdout : stderr);
    static const char str[] = {
            "SYNOPSIS\n"
            "  IUPAC FASTA extraction from BAM files to stdout\n"
            "USAGE\n"
            "  %s [options] --b 1.bam --n ind1 --b 2.bam --n ind2 [...]\n"
            "  OR\n"
            "  %s [options] --bamList samples.txt\n"
            "OPTIONS\n"
            "  --b\t\tBAM files\n"
            "  --n\t\tSequence IDs\n"
            "  --bamList\tFile containing BAM paths and Sample IDs (one per line: <BAM> <SAMPLE>)\n"
            "  --r\t\tRegion ('chr:start-end' or 'chr start end') coordinates are 1-based\n"
            "  --regionList\tFile containing regions (one per line: 'chr:start-end' or 'chr start end')\n"
            "  --minMQ\tMinimum mapping quality (default: 0)\n"
            "  --minBQ\tMinimum base quality (default: 0)\n"
            "  --minC\tMinimum coverage (default: 0)\n"
            "  --maxC\tMaximum coverage (default: 9999)\n"
            "  --iupacRatio\tIUPAC ratio (default: 0.25)\n"
            "  --incMQ\tInclude missing mapping quality value 255 (default: False)\n"
            "  --incBQ\tInclude missing base quality value 255 (default: False)\n"
            "  --help\tShow this help\n"
            "  --version\tPrint version and exit\n"
            "  --debug\tDebug\n"
            "\n"
            "NOTE: Multiple BAM files (--b/--bamList) and regions (--r/--regionList)\n"
            "      are processed in the order provided.\n"
            "      Multiple regions are concatenated into a single output sequence,\n"
            "      enabling direct extraction of combined intervals such as exons from a GTF annotation.\n"
            "\n"
            "EXENAME\n"
            "  %s\n"
            "VERSION\n"
            "  %s\n"
            "URL\n"
            "  %s\n"};
    (void)fprintf(out, str, program_name, program_name, EXENAME, VERSION, GITHUB_URL);
    exit(retcode);
}

void show_version(
    ) {
    
    std::cout << EXENAME << std::endl;
    std::cout << "Version: " << VERSION << std::endl;
    exit(EXIT_SUCCESS);
}

int main(
    int argc,
    char** argv) {

    std::vector<std::string> bamFiles;
    std::vector<std::string> sequenceIds;
    std::vector<std::string> regions;
    std::string bamList = "";
    std::string regionList = "";
    int minMapQuality = 0;
    int maxMapQuality = 254;
    int minBaseQuality = 0;
    int maxBaseQuality = 254;
    int minCoverage = 0;
    int maxCoverage = 9999;
    double iupacRatio = 0.25;
    int debug = 0;

    static struct option long_options[] = {
            {"b", required_argument, NULL, 0 },
            {"n", required_argument, NULL, 1 },
            {"r", required_argument, NULL, 2 },
            {"minMQ", required_argument, NULL, 3 },
            {"minBQ", required_argument, NULL, 4 },
            {"minC", required_argument, NULL, 5 },
            {"maxC", required_argument, NULL, 6 },
            {"iupacRatio", required_argument, NULL, 7 },
            {"incMQ", no_argument, NULL, 8 },
            {"incBQ", no_argument, NULL, 9 },
            {"help", no_argument, NULL, 10 },
            {"version", no_argument, NULL, 11 },
            {"debug", no_argument, NULL, 12 },
            {"bamList", required_argument, NULL, 13 },
            {"regionList", required_argument, NULL, 14 },
            {NULL, 0, NULL, 0 }
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                if (strcmp(long_options[option_index].name, "b") == 0) {
                    bamFiles.push_back(optarg);
                }
                break;
            case 1:
                if (strcmp(long_options[option_index].name, "n") == 0) {
                    sequenceIds.push_back(optarg);
                }
                break;
            case 2:
                if (strcmp(long_options[option_index].name, "r") == 0) {
                    regions.push_back(optarg);
                }
                break;
            case 3:
                if (strcmp(long_options[option_index].name, "minMQ") == 0) {
                    minMapQuality = std::stoi(optarg);
                }
                break;
            case 4:
                if (strcmp(long_options[option_index].name, "minBQ") == 0) {
                    minBaseQuality = std::stoi(optarg);
                }
                break;
            case 5:
                if (strcmp(long_options[option_index].name, "minC") == 0) {
                    minCoverage = std::stoi(optarg);
                }
                break;
            case 6:
                if (strcmp(long_options[option_index].name, "maxC") == 0) {
                    maxCoverage = std::stoi(optarg);
                }
                break;
            case 7:
                if (strcmp(long_options[option_index].name, "iupacRatio") == 0) {
                    iupacRatio = std::stod(optarg);
                }
                break;
            case 8:
                if (strcmp(long_options[option_index].name, "incMQ") == 0) {
                    maxMapQuality = 255;
                }
                break;
            case 9:
                if (strcmp(long_options[option_index].name, "incBQ") == 0) {
                    maxBaseQuality = 255;
                }
                break;
            case 10:
                if (strcmp(long_options[option_index].name, "help") == 0) {
                    show_help(argv[0], EXIT_SUCCESS);
                }
                break;
            case 11:
                if (strcmp(long_options[option_index].name, "version") == 0) {
                    show_version();
                }
                break;
            case 12:
                if (strcmp(long_options[option_index].name, "debug") == 0) {
                    debug = 1;
                }
                break;
            case 13:
                if (strcmp(long_options[option_index].name, "bamList") == 0) {
                    bamList = optarg;
                }
                break;
            case 14:
                if (strcmp(long_options[option_index].name, "regionList") == 0) {
                    regionList = optarg;
                }
                break;
            case '?':
                std::cerr << "Invalid option or missing argument" << std::endl;
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
                return 1;
        }
    }

    // Check if required options are provided
    if (!bamList.empty()) {
        std::ifstream infile(bamList);
        if (!infile.is_open()) {
            std::cerr << "Error: Could not open bam list file " << bamList << std::endl;
            return 1;
        }
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty()) {
                continue;
            }
            if (line[0] == '#') {
                continue;
            }
            std::stringstream ss(line);
            std::string bam, sample;
            if (ss >> bam >> sample) {
                bamFiles.push_back(bam);
                sequenceIds.push_back(sample);
            } else {
                std::cerr << "Warning: Skipping malformed line in bam list file: " << line << std::endl;
            }
        }
        infile.close();
    }
    if (!regionList.empty()) {
        std::ifstream infile(regionList);
        if (!infile.is_open()) {
            std::cerr << "Error: Could not open region list file " << regionList << std::endl;
            return 1;
        }
        std::string line;
        while (std::getline(infile, line)) {
            if (line.empty()) {
                continue;
            }
            if (line[0] == '#') {
                continue;
            }
            regions.push_back(line);
        }
        infile.close();
    }
    if (bamFiles.empty()) {
        std::cerr << "Error: At least one BAM file must be specified with the --b option (see --help)." << std::endl;
        return 1;
    }
    if (regions.empty()) {
        std::cerr << "Error: At least one region must be specified with the --r option (see --help)." << std::endl;
        return 1;
    }
    if (!sequenceIds.empty() && sequenceIds.size() != bamFiles.size()) {
        std::cerr << "Error: Number of sequence IDs (--n) must match number of BAM files (--b) (see --help)." << std::endl;
        return 1;
    }

    // Parse region(s) string to extract chromosome, start position, and end position
    std::vector<std::string> chromosomes;
    std::vector<int> startPositions;
    std::vector<int> endPositions;
    std::tuple<std::string, int, int> regionSplit;
    std::string chromosome = "";
    int startPos = 0;
    int endPos = 0;
    regionSplit = std::make_tuple(chromosome, startPos, endPos);
    for (const std::string& region : regions) {
        if (debug == 1) {
            std::cout << "region: " << region << std::endl;
        }
        try {
            regionSplit = parseRegion(region);
        } catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
            // Exit
            return EXIT_FAILURE;
        }

        if (debug == 1) {
            std::cout << "chromosome: " << std::get<0>(regionSplit) << std::endl;
            std::cout << "startPos: " << std::get<1>(regionSplit) << std::endl;
            std::cout << "endPos: " << std::get<2>(regionSplit) << std::endl;
        }
        chromosomes.push_back(std::get<0>(regionSplit));
        startPositions.push_back(std::get<1>(regionSplit));
        endPositions.push_back(std::get<2>(regionSplit));
    }

    // Use BAM file names as sequence IDs if none are provided
    if (sequenceIds.empty()) {
        sequenceIds = bamFiles;
    }
    if (debug == 1) {
        for (const auto &sequenceId : sequenceIds) {
            std::cout << "seqID: " <<  sequenceId << std::endl;
        }
    }

    // Vector to store results for each BAM file
    std::vector<std::vector<std::vector<std::vector<int>>>> allCounts;
    allCounts.resize(regions.size());
    std::vector<std::vector<std::vector<std::string>>> allIupacStrings;
    allIupacStrings.resize(regions.size());

    for (size_t r = static_cast<size_t>(0); r < regions.size(); ++r) {
        bool regionResolved = false;
        for (const auto &bamFile : bamFiles) {
            BamExtractionResult bamresult;
            try {
                bamresult = extractCountsFromBam(
                    bamFile,
                    chromosomes[r],
                    startPositions[r],
                    endPositions[r],
                    minMapQuality,
                    maxMapQuality,
                    minBaseQuality,
                    maxBaseQuality,
                    minCoverage,
                    maxCoverage,
                    debug);
            } catch (const std::runtime_error& e) {
                std::cerr << e.what() << std::endl;
                // Exit
                return EXIT_FAILURE;
            }
            if (!regionResolved) {
                startPositions[r] = bamresult.startPos;
                endPositions[r] = bamresult.endPos;
                regionResolved = true;
            }

            std::vector<std::vector<int>> counts = bamresult.counts;

            // Vector to store merged IUPAC characters for current BAM file
            std::vector<std::string> iupacStrings;
            for (const auto &count: counts) {
                char iupacChar = getIupac(count, iupacRatio); // Pass appropriate arguments
                iupacStrings.push_back(std::string(static_cast<size_t>(1), iupacChar));
            }
            allCounts[r].push_back(counts);
            allIupacStrings[r].push_back(iupacStrings);
        }
    }
    if (debug == 1) {
        printAllCounts(bamFiles, chromosomes, startPositions, endPositions, allCounts);
    }
    // Generate FASTA-formatted sequences
    std::string fastaSequences = generateFasta(allIupacStrings, sequenceIds, chromosomes, startPositions, endPositions);

    // Output FASTA-formatted sequences to stdout
    std::cout << fastaSequences;

    return EXIT_SUCCESS;
}
