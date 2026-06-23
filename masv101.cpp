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

// --- Global Structures & Config ---
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

// --- Helper Functions ---
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

std::vector<SequenceData> dereplicate_fasta(const std::string& filename) {
    std::unordered_map<std::string, int> seq_counter;
    std::vector<std::string> order;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << "\n";
        exit(1);
    }

    std::string line, current_seq = "";
    int total_count = 0;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                if (seq_counter.find(current_seq) == seq_counter.end()) order.push_back(current_seq);
                seq_counter[current_seq]++;
                total_count++;
                current_seq.clear();
            }
        } else {
            for (char &c : line) c = toupper(c);
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        if (seq_counter.find(current_seq) == seq_counter.end()) order.push_back(current_seq);
        seq_counter[current_seq]++;
        total_count++;
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
        char title_buf[64];
        sprintf(title_buf, "Uniq%04zu;size=%d;", i + 1, sorted_seqs[i].second);
        
        SequenceData sd;
        sd.title = std::string(title_buf);
        sd.seq = sorted_seqs[i].first;
        sd.size = sorted_seqs[i].second;
        sd.length = static_cast<int>(sd.seq.length());
        data.push_back(std::move(sd));
    }

    std::cout << "Total input sequences: " << total_count << "\n";
    std::cout << "Dereplicated sequences: " << data.size() << "\n";
    return data;
}

// Global outputs protected by a single Mutex (Only used when a match is actually found)
std::mutex output_mutex;
std::unordered_map<std::string, std::vector<std::string>> filtered_map;
std::unordered_map<std::string, std::string> variants_map;
std::unordered_map<std::string, std::string> noise_map;
std::unordered_map<std::string, int> weights_map;

struct NoiseMeta {
    int target_idx;
    int p_sum;
    int im_sum;
    int len_diff;
};

int main(int argc, char* argv[]) {
    std::string input_fasta = "";
    int fx = 1;
    float ax = 1.45f;
    bool save_spurious = false;
    int num_threads = 2;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) input_fasta = argv[++i];
        else if (arg == "-f" && i + 1 < argc) fx = std::stoi(argv[++i]);
        else if (arg == "-a" && i + 1 < argc) ax = std::stof(argv[++i]);
        else if (arg == "-s" && i + 1 < argc) save_spurious = (std::string(argv[++i]) == "True");
        else if (arg == "-t" && i + 1 < argc) num_threads = std::stoi(argv[++i]);
    }

    if (input_fasta.empty()) {
        std::cerr << "Usage: " << argv[0] << " -i <input.fasta> [-f freedom] [-a abundance] [-s True/False] [-t threads]\n";
        return 1;
    }

    std::cout << "\nMASV v.1.0.1 (Multi-Threaded C++ Engine)\n";
    auto start_time = std::chrono::high_resolution_clock::now();

    auto dataset = dereplicate_fasta(input_fasta);
    size_t N = dataset.size();
    if (N == 0) return 0;

    std::cout << "Generating K-mer feature profiles across " << num_threads << " thread(s)...\n";
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

    std::cout << "Building algorithmic lookup buckets...\n";
    std::unordered_map<int, std::vector<int>> length_buckets;
    for (int i = 0; i < static_cast<int>(N); ++i) {
        length_buckets[dataset[i].length].push_back(i);
    }

    std::cout << "Identifying sequence variants...\n";
    
    // THE FIX: Lock-Free Atomic Array
    std::unique_ptr<std::atomic<bool>[]> active(new std::atomic<bool>[N]);
    for(size_t i = 0; i < N; ++i) {
        active[i].store(true, std::memory_order_relaxed);
    }

    std::atomic<size_t> main_index(0);
    std::atomic<size_t> progress_counter(0);
    std::vector<std::thread> worker_threads;

    for (int t = 0; t < num_threads; ++t) {
        worker_threads.emplace_back([&]() {
            
            // THE FIX: Pre-allocate memory ONCE per thread, outside the loop
            std::vector<int> candidates;
            candidates.reserve(10000); 
            std::vector<size_t> flagged_noise_indices;
            flagged_noise_indices.reserve(10000);
            std::vector<NoiseMeta> noise_metadata;
            noise_metadata.reserve(10000);

            while (true) {
                size_t i = main_index.fetch_add(1, std::memory_order_relaxed);
                if (i >= N) break;

                // Ultra-fast lock-free read
                if (!active[i].load(std::memory_order_relaxed)) {
                    progress_counter.fetch_add(1, std::memory_order_relaxed);
                    continue;
                }

                candidates.clear();
                flagged_noise_indices.clear();
                noise_metadata.clear();

                int L_i = dataset[i].length;
                for (int target_L = L_i - 1; target_L <= L_i + 1; ++target_L) {
                    auto it = length_buckets.find(target_L);
                    if (it != length_buckets.end()) {
                        for (int idx : it->second) {
                            if (idx > static_cast<int>(i)) candidates.push_back(idx);
                        }
                    }
                }

                int loop_count = 0;

                for (int target_idx : candidates) {
                    // Ultra-fast lock-free read
                    if (!active[target_idx].load(std::memory_order_relaxed)) continue;

                    if ((static_cast<float>(dataset[i].size) / dataset[target_idx].size) < ax) continue;

                    int p_sum = 0;
                    int im_sum = 0;
                    bool failed = false;

                    for (int k = 0; k < 16; ++k) {
                        int diff = std::abs(dataset[i].kmers[k] - dataset[target_idx].kmers[k]);
                        if (k == 0 || k == 5 || k == 10 || k == 15) {
                            p_sum += diff;
                            if (p_sum > 4 * fx) { failed = true; break; } 
                        } else {
                            im_sum += diff;
                            if (im_sum > 4 * fx) { failed = true; break; } 
                        }
                    }
                    if (failed || (im_sum - p_sum) > 2 * fx) continue;

                    flagged_noise_indices.push_back(target_idx);
                    int len_diff = std::abs(dataset[i].length - dataset[target_idx].length);
                    noise_metadata.push_back({target_idx, p_sum, im_sum, len_diff});
                    loop_count++;
                }

                // Ultra-fast lock-free write (Deactivate noise sequences immediately)
                for (size_t idx : flagged_noise_indices) {
                    active[idx].store(false, std::memory_order_relaxed);
                }

                // Only lock when writing confirmed variants to the global dictionary
                if (!noise_metadata.empty() || dataset[i].size >= 2) {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    weights_map[dataset[i].title] = loop_count;

                    for (const auto& meta : noise_metadata) {
                        std::string t_title = dataset[meta.target_idx].title;
                        noise_map[t_title] = dataset[meta.target_idx].seq;
                        filtered_map[t_title] = {t_title, dataset[i].title, "NOISY VARIANT", std::to_string(meta.p_sum), std::to_string(meta.im_sum), std::to_string(meta.len_diff)};
                    }

                    if (dataset[i].size >= 2) {
                        variants_map[dataset[i].title] = dataset[i].seq;
                        filtered_map[dataset[i].title] = {dataset[i].title, "*", "VARIANT", "*", "*", "*"};
                    }
                }

                size_t done = progress_counter.fetch_add(1, std::memory_order_relaxed) + 1;
                if (done % std::max((size_t)1, N / 100) == 0 || done == N) {
                    std::cout << "\rProcessing updates: " << (done * 100 / N) << "% Completed." << std::flush;
                }
            }
        });
    }
    for (auto& th : worker_threads) th.join();

    for (size_t i = 0; i < N; ++i) {
        if (filtered_map.find(dataset[i].title) == filtered_map.end()) {
            if (save_spurious) {
                variants_map[dataset[i].title] = dataset[i].seq;
            } else {
                noise_map[dataset[i].title] = dataset[i].seq;
            }
            filtered_map[dataset[i].title] = {dataset[i].title, "*", "SPURIOUS VARIANT", "*", "*", "*"};
        }
    }

    std::cout << "\nWriting data outputs to files...\n";
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
            var_file << ">" << it->first << "neighbors=" << weights_map[it->first] << ";\n" << it->second << "\n";
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
    
    std::cout << "\nDetected sequence variants: " << variants_map.size() << "\n";
    std::cout << "Script execution time: " << total_time.count() << " seconds.\n\n";

    return 0;
}