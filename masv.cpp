#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <cmath>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <memory>
#include <iomanip>

const int PERFECT_IDX[] = {0, 5, 10, 15};
const int IMPERFECT_IDX[] = {1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14};

inline int base_to_val(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

inline int kmer_to_idx(char b1, char b2) {
    int v1 = base_to_val(b1);
    int v2 = base_to_val(b2);
    if (v1 == -1 || v2 == -1) return -1;
    return (v1 << 2) | v2;
}

struct SequenceData {
    std::string title;
    std::string seq;
    int size;
    int length;
    std::vector<int> kmers;
};

struct AdjEdge {
    int target;
    int p_sum;
    int im_sum;
    int len_diff;
};

std::vector<int> generate_kmer_vector(const std::string& seq, const std::string& title) {
    std::vector<int> vector(16, 0);
    for (char c : seq) {
        int val = base_to_val(c);
        if (val == -1) {
            std::cerr << "\nError processing " << title << ": Ambiguous nucleotide detected.\n";
            exit(1);
        }
        vector[(val << 2) | val]++; 
    }
    for (size_t i = 0; i < seq.length() - 1; ++i) {
        int idx = kmer_to_idx(seq[i], seq[i + 1]);
        if (idx != -1) vector[idx]++;
    }
    return vector;
}

// Count one read per record, or the header's size= value when it is present.
std::vector<SequenceData> parse_fasta(const std::string& filename) {
    std::unordered_map<std::string, int> seq_counter;
    std::unordered_map<std::string, std::string> original_titles;
    std::vector<std::string> order;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        exit(1);
    }

    std::string line, current_seq = "", current_title = "";
    int total_count = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                if (seq_counter.find(current_seq) == seq_counter.end()) {
                    order.push_back(current_seq);
                    original_titles[current_seq] = current_title;
                }
                
                size_t size_pos = current_title.find("size=");
                if (size_pos != std::string::npos) {
                    size_t end_pos = current_title.find(';', size_pos);
                    if (end_pos == std::string::npos) end_pos = current_title.length();
                    int parsed_size = std::stoi(current_title.substr(size_pos + 5, end_pos - (size_pos + 5)));
                    seq_counter[current_seq] += parsed_size;
                    total_count += parsed_size;
                } else {
                    seq_counter[current_seq]++;
                    total_count++;
                }
                current_seq.clear();
            }
            current_title = line.substr(1);
        } else {
            for (char &c : line) c = toupper(c);
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        if (seq_counter.find(current_seq) == seq_counter.end()) {
            order.push_back(current_seq);
            original_titles[current_seq] = current_title;
        }
        size_t size_pos = current_title.find("size=");
        if (size_pos != std::string::npos) {
            size_t end_pos = current_title.find(';', size_pos);
            if (end_pos == std::string::npos) end_pos = current_title.length();
            int parsed_size = std::stoi(current_title.substr(size_pos + 5, end_pos - (size_pos + 5)));
            seq_counter[current_seq] += parsed_size;
            total_count += parsed_size;
        } else {
            seq_counter[current_seq]++;
            total_count++;
        }
    }
    file.close();

    std::vector<std::pair<std::string, int>> sorted_seqs;
    sorted_seqs.reserve(order.size());
    for (const auto& seq : order) {
        sorted_seqs.push_back({seq, seq_counter[seq]});
    }
    
    std::sort(sorted_seqs.begin(), sorted_seqs.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    std::vector<SequenceData> data;
    data.reserve(sorted_seqs.size());
    for (size_t i = 0; i < sorted_seqs.size(); ++i) {
        SequenceData sd;
        std::string orig_title = original_titles[sorted_seqs[i].first];
        
        if (orig_title.find("size=") != std::string::npos) {
            sd.title = orig_title; 
        } else {
            sd.title = orig_title + ";size=" + std::to_string(sorted_seqs[i].second) + ";";
        }
        
        sd.seq = sorted_seqs[i].first;
        sd.size = sorted_seqs[i].second;
        sd.length = static_cast<int>(sd.seq.length());
        data.push_back(std::move(sd));
    }

    std::cout << total_count << " seqs, " << data.size() << " uniques\n";
    return data;
}

// Classification results are keyed by sequence title.
std::unordered_map<std::string, std::vector<std::string>> filtered_map;
std::unordered_map<std::string, std::string> variants_map;
std::unordered_map<std::string, std::string> noise_map;
std::unordered_map<std::string, int> weights_map;

int main(int argc, char* argv[]) {
    std::string input_fasta = "";
    float ax = 2.0f;
    bool save_spurious = false;
    int num_threads = 2;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) input_fasta = argv[++i];
        else if (arg == "-a" && i + 1 < argc) ax = std::stof(argv[++i]);
        else if (arg == "-s" && i + 1 < argc) save_spurious = (std::string(argv[++i]) == "True");
        else if (arg == "-t" && i + 1 < argc) num_threads = std::stoi(argv[++i]);
    }

    if (input_fasta.empty()) {
        std::cerr << "Usage: " << argv[0] << " -i <input.fasta> [-a abundance] [-s True/False] [-t threads]\n";
        return 1;
    }

    std::cout << "\nMASV 2.0.0\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    auto dataset = parse_fasta(input_fasta);
    size_t N = dataset.size();
    if (N == 0) return 0;

    std::cout << "Counting k-mers...\n";
    std::vector<std::thread> kmer_threads;
    std::atomic<size_t> kmer_index(0);

    for (int t = 0; t < num_threads; ++t) {
        kmer_threads.emplace_back([&]() {
            while (true) {
                size_t curr = kmer_index.fetch_add(1, std::memory_order_relaxed);
                if (curr >= N) break;
                dataset[curr].kmers = generate_kmer_vector(dataset[curr].seq, dataset[curr].title);
            }
        });
    }
    for (auto& th : kmer_threads) th.join();

    std::cout << "Indexing...\n";
    std::unordered_map<int, std::vector<int>> length_buckets;
    for (int i = 0; i < static_cast<int>(N); ++i) {
        length_buckets[dataset[i].length].push_back(i);
    }

    std::cout << "Denoising\n";
    std::vector<int> parent(N, -1);
    std::vector<AdjEdge> best_edge(N);

    std::atomic<size_t> main_index(0);
    std::atomic<size_t> progress_counter(0);
    std::vector<std::thread> worker_threads;

    for (int t = 0; t < num_threads; ++t) {
        worker_threads.emplace_back([&]() {
            std::vector<int> candidates;
            candidates.reserve(10000);

            while (true) {
                size_t i = main_index.fetch_add(1, std::memory_order_relaxed);
                if (i >= N) break;

                candidates.clear();
                int L_i = dataset[i].length;

                // Candidates: index j < i (at least as abundant) and length within +/-1 nt.
                for (int target_L = L_i - 1; target_L <= L_i + 1; ++target_L) {
                    auto it = length_buckets.find(target_L);
                    if (it != length_buckets.end()) {
                        for (int j : it->second) {
                            if (j < static_cast<int>(i)) {
                                candidates.push_back(j);
                            }
                        }
                    }
                }

                std::sort(candidates.begin(), candidates.end());

                int best_j = -1;
                AdjEdge temp_edge = {-1, 0, 0, 0};

                for (int j : candidates) {
                    if ((static_cast<float>(dataset[j].size) / dataset[i].size) < ax) continue;

                    int p_sum = 0, im_sum = 0;
                    bool failed = false;

                    for (int k = 0; k < 16; ++k) {
                        int diff = std::abs(dataset[i].kmers[k] - dataset[j].kmers[k]);
                        if (k == 0 || k == 5 || k == 10 || k == 15) {
                            p_sum += diff;
                            if (p_sum > 4) { failed = true; break; } 
                        } else {
                            im_sum += diff;
                            if (im_sum > 4) { failed = true; break; } 
                        }
                    }
                    
                    if (!failed && (im_sum - p_sum) <= 2) {
                        best_j = j;
                        temp_edge = {j, p_sum, im_sum, std::abs(L_i - dataset[j].length)};
                        break; 
                    }
                }

                if (best_j != -1) {
                    parent[i] = best_j;
                    best_edge[i] = temp_edge;
                }

                size_t done = progress_counter.fetch_add(1, std::memory_order_relaxed) + 1;
                if (done % std::max((size_t)1, N / 100) == 0 || done == N) {
                    float pct = (static_cast<float>(done) * 100.0f) / N;
                    std::cout << "\r  " << std::fixed << std::setprecision(1) << pct << "% " << std::flush;
                }
            }
        });
    }
    for (auto& th : worker_threads) th.join();
    std::cout << "\n";

    std::vector<bool> is_variant(N, true);
    std::vector<int> direct_noise_mass(N, 0);

    for (size_t i = 0; i < N; ++i) {
        if (parent[i] != -1) {
            is_variant[i] = false;
            direct_noise_mass[parent[i]] += dataset[i].size;
        }
    }

    std::cout << "Writing outputs\n";
    for (size_t i = 0; i < N; ++i) {
        std::string t_title = dataset[i].title;

        if (is_variant[i]) {
            if (!save_spurious && dataset[i].size == 1 && direct_noise_mass[i] == 0) {
                noise_map[t_title] = dataset[i].seq;
                filtered_map[t_title] = {t_title, "*", "SPURIOUS VARIANT", "*", "*", "*"};
            } else {
                variants_map[t_title] = dataset[i].seq;
                weights_map[t_title] = direct_noise_mass[i]; 
                filtered_map[t_title] = {t_title, "*", "VARIANT", "*", "*", "*"};
            }
        } else {
            noise_map[t_title] = dataset[i].seq;
            int p_idx = parent[i];
            std::string rep_title = dataset[p_idx].title;
            
            filtered_map[t_title] = {
                t_title, 
                rep_title, 
                "NOISY VARIANT", 
                std::to_string(best_edge[i].p_sum), 
                std::to_string(best_edge[i].im_sum), 
                std::to_string(best_edge[i].len_diff)
            };
        }
    }

    std::ofstream tab_file("asv_tab.txt");
    tab_file << "title\tclosest neighbor(s)\tdescription\tperfect k-mer\timperfect k-mer\tlength difference\n";
    for (size_t i = 0; i < N; ++i) {
        auto it = filtered_map.find(dataset[i].title);
        if (it != filtered_map.end()) {
            tab_file << it->second[0] << "\t" << it->second[1] << "\t" << it->second[2] << "\t" 
                     << it->second[3] << "\t" << it->second[4] << "\t" << it->second[5] << "\n";
        }
    }
    tab_file.close();

    std::ofstream var_file("variants.fa");
    for (size_t i = 0; i < N; ++i) {
        auto it = variants_map.find(dataset[i].title);
        if (it != variants_map.end()) {
            var_file << ">" << it->first << "noise=" << weights_map[it->first] << ";\n" << it->second << "\n";
        }
    }
    var_file.close();

    std::ofstream noise_file("noise.fa");
    for (size_t i = 0; i < N; ++i) {
        auto it = noise_map.find(dataset[i].title);
        if (it != noise_map.end()) {
            noise_file << ">" << it->first << "\n" << it->second << "\n";
        }
    }
    noise_file.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_time = end_time - start_time;
    
    std::cout << variants_map.size() << " ASVs\n";
    std::cout << "Time " << std::fixed << std::setprecision(1) << total_time.count() << "s\n\n";

    return 0;
}
