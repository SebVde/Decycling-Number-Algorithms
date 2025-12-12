import datetime
import multiprocessing

import networkx as nx
from pyvis.network import Network
from exact_mif import get_decycling_number_mif
from naive import get_decycling_number
from exact_mif_v2 import get_decycling_number_mif_v2
from exact_mif_v3 import get_decycling_number_mif_v3
from bafna_fvs import get_decycling_number_2_approx
from constr_heur import get_decycling_number_constr_heur
import re
import io
import numpy as np
import time
import os
from queue import Empty as QueueEmpty

all_matrices = []


def parse_adj_matrices(file):
    graphs = []
    try:
        with open(file, "r") as f:
            content = f.read()

        matrix_blocks = re.split(r"\n\s*\n", content.strip())

        for i, block in enumerate(matrix_blocks):
            if not block.strip():
                continue

            try:
                all_matrices.append(block)
                np_matrix = np.loadtxt(io.StringIO(block), dtype=int)
                G = nx.from_numpy_array(np_matrix)
                graphs.append(G)
                if len(graphs) == 20:
                    break

            except Exception as e:
                print(f"Error at block {i + 1}: {e}")

    except FileNotFoundError:
        print(f"Error : The file '{file}' was not found.")
        return []

    return graphs


def run_method_in_process(method_func, graph_data, result_queue):
    try:
        start_time = time.perf_counter()
        method_func(graph_data)
        end_time = time.perf_counter()
        result_queue.put(end_time - start_time)
    except Exception as e:
        result_queue.put(e)


def benchmark_graph_methods(
    directory_path,
    methods_list,
    timeout_seconds,
    output_filename="benchmark_results.txt",
):
    benchmark_start_time = time.perf_counter()
    file_pattern = re.compile(r"v_(\d+)_D_(\d+)\.mat")

    if not os.path.isdir(directory_path):
        print(f"Error : Folder '{directory_path}' does not exist.")
        return

    # files_to_process = sorted(
    #     f for f in os.listdir(directory_path) if file_pattern.match(f)
    # )

    files_to_process = [f for f in os.listdir(directory_path)]

    if not files_to_process:
        print(f"No file to process has been found in folder '{directory_path}'.")
        return

    try:
        with open(output_filename, "w", encoding="utf-8") as f_out:

            f_out.write("=" * 80 + "\n")
            f_out.write(f"BENCHMARK RESULTS\n")
            f_out.write(f"Date: {time.ctime()}\n")
            f_out.write(f"Treated folder: {os.path.abspath(directory_path)}\n")
            f_out.write(f"Timeout per graph: {timeout_seconds}s\n")
            f_out.write("=" * 80 + "\n\n")

            for filename in files_to_process:

                file_header_line = f"Treating: {filename}"
                file_separator = "=" * 80

                print(f"\n{file_separator}")
                print(file_header_line)
                print(f"{file_separator}")

                f_out.write(f"{file_separator}\n")
                f_out.write(f"{file_header_line}\n")
                f_out.write(f"{file_separator}\n\n")

                filepath = os.path.join(directory_path, filename)

                try:
                    graphs_in_file = parse_adj_matrices(filepath)
                    total_graphs = len(graphs_in_file)
                    if total_graphs == 0:
                        print(f"No graphs loaded from '{filename}'. Skipping.")
                        continue
                    msg = f"  File loaded: {total_graphs} graphs."
                    print(msg)
                    f_out.write(msg + "\n")

                except Exception as e:
                    print(f"Error loading graphs from '{filename}': {e}")
                    continue

                results = {
                    m.__name__: {"times": [], "timeouts": 0, "errors": 0}
                    for m in methods_list
                }

                for i, graph in enumerate(graphs_in_file):
                    print(f"  Test... Graph {i + 1}/{total_graphs}", end="\r")

                    for method in methods_list:
                        method_name = method.__name__

                        run_times = []
                        timeout_occurred = False
                        error_occurred = False

                        for _ in range(3):
                            result_queue = multiprocessing.Queue()

                            p = multiprocessing.Process(
                                target=run_method_in_process,
                                args=(method, graph, result_queue),
                            )
                            p.start()
                            p.join(timeout=timeout_seconds)

                            if p.is_alive():
                                p.terminate()
                                p.join()
                                timeout_occurred = True
                                break
                            elif p.exitcode != 0:
                                error_occurred = True
                                break
                            else:
                                try:
                                    result = result_queue.get_nowait()
                                    if isinstance(result, Exception):
                                        error_occurred = True
                                        break
                                    else:
                                        run_times.append(result)
                                except QueueEmpty:
                                    error_occurred = True
                                    break

                        if timeout_occurred:
                            results[method_name]["timeouts"] += 1
                        elif error_occurred:
                            results[method_name]["errors"] += 1
                        elif run_times:
                            median_time = np.median(run_times)
                            results[method_name]["times"].append(median_time)

                print(f"\n  Finished tests for {filename}.\n")
                f_out.write("\n  Finished tests. Results incoming...\n\n")

                col_methode = max(len(m.__name__) for m in methods_list) + 2
                col_methode = max(col_methode, 10)  # min 10
                col_time = 12
                col_pct = 12

                table_header = (
                    f"| {'Method':<{col_methode}} "
                    f"| {'Min (s)':>{col_time}} "
                    f"| {'Max (s)':>{col_time}} "
                    f"| {'Mean (s)':>{col_time}} "
                    f"| {'Median (s)':>{col_time}} "
                    f"| {'Success (%)':>{col_pct}} "
                    f"| {'Timeout (%)':>{col_pct}} "
                    f"| {'Errors (%)':>{col_pct}} |"
                )
                table_separator = (
                    f"|{'-' * (col_methode + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_pct + 1)}|"
                )

                print(table_header)
                print(table_separator)
                f_out.write(table_header + "\n")
                f_out.write(table_separator + "\n")

                for method_name, data in results.items():
                    times_sec = np.array(data["times"])
                    num_success = len(times_sec)
                    num_timeouts = data["timeouts"]
                    num_errors = data["errors"]

                    timeout_percent = (num_timeouts / total_graphs) * 100
                    error_percent = (num_errors / total_graphs) * 100
                    success_percent = (num_success / total_graphs) * 100

                    if num_success > 0:
                        s_min = f"{np.min(times_sec):.6f}"
                        s_max = f"{np.max(times_sec):.6f}"
                        s_mean = f"{np.mean(times_sec):.6f}"
                        s_median = f"{np.median(times_sec):.6f}"
                    else:
                        s_min = s_max = s_mean = s_median = "---"

                    row_data = (
                        f"| {method_name:<{col_methode}} "
                        f"| {s_min:>{col_time}} "
                        f"| {s_max:>{col_time}} "
                        f"| {s_mean:>{col_time}} "
                        f"| {s_median:>{col_time}} "
                        f"| {success_percent:>{col_pct - 3}.2f} % "
                        f"| {timeout_percent:>{col_pct - 3}.2f} % "
                        f"| {error_percent:>{col_pct - 3}.2f} % |"
                    )

                    print(row_data)
                    f_out.write(row_data + "\n")

                f_out.write("\n")

            benchmark_end_time = time.perf_counter()
            total_seconds = benchmark_end_time - benchmark_start_time
            formatted_time = str(datetime.timedelta(seconds=total_seconds))
            print(f"\n{file_separator}")
            print(f"Benchmark over. Results written in '{output_filename}'.")
            print(f"Total execution time: {formatted_time}")
            print(f"{file_separator}")
            f_out.write(f"\n{file_separator}\n")
            f_out.write("Benchmark over.\n")
            f_out.write(f"Total execution time: {formatted_time}\n")

    except IOError as e:
        print(f"Error: Impossible to write in '{output_filename}': {e}")

    except Exception as e:
        print(f"Error : {e}")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    TIMEOUT_MAX = 600

    METHODS = [
        get_decycling_number_mif,
        get_decycling_number_mif_v2,
        get_decycling_number_mif_v3,
    ]

    # benchmark_graph_methods(
    #     directory_path="Benchmark graphs/diameter",
    #     methods_list=METHODS,
    #     timeout_seconds=TIMEOUT_MAX,
    #     output_filename="benchmark_results_diameter.txt",
    # )
    #
    # benchmark_graph_methods(
    #     directory_path="Benchmark graphs/chromatic number",
    #     methods_list=METHODS,
    #     timeout_seconds=TIMEOUT_MAX,
    #     output_filename="benchmark_results_chrom_number.txt",
    # )
    #
    # benchmark_graph_methods(
    #     directory_path="Benchmark graphs/density",
    #     methods_list=METHODS,
    #     timeout_seconds=TIMEOUT_MAX,
    #     output_filename="benchmark_results_density.txt",
    # )

    benchmark_graph_methods(
        directory_path="Benchmark graphs/domination number",
        methods_list=METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="benchmark_results_domination_number.txt",
    )

    benchmark_graph_methods(
        directory_path="Benchmark graphs/girth",
        methods_list=METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="benchmark_results_girth.txt",
    )

    benchmark_graph_methods(
        directory_path="Benchmark graphs/longest induced cycle",
        methods_list=METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="benchmark_results_long_ind_cyc.txt",
    )

    benchmark_graph_methods(
        directory_path="Benchmark graphs/radius",
        methods_list=METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="benchmark_results_radius.txt",
    )

    benchmark_graph_methods(
        directory_path="Benchmark graphs/treewidth",
        methods_list=METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="benchmark_results_treewidth.txt",
    )

    benchmark_graph_methods(
        directory_path="Benchmark graphs/vertex connectivity",
        methods_list=METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="benchmark_results_vertex_con.txt",
    )

# graphs = parse_adj_matrices("Benchmark graphs/equal_to_1.mat")
# for graph in graphs:
#     print("Decycling number (naive):", get_decycling_number(graph))
#     print("Decycling number (MIF):", get_decycling_number_mif(graph))
#     print("Decycling number (MIF v2):", get_decycling_number_mif_v2(graph))
#     print("Decycling number (MIF v3):", get_decycling_number_mif_v3(graph))
#     print("-----")

# nt = nx.erdos_renyi_graph(900, 0.9)
# Densité faible -> plus lent
# print(get_decycling_number_constr_heur(nt))

# nt = nx.Graph()
# nt.add_nodes_from(["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8"])
# nt.add_edges_from(
#     [
#         ("V1", "V2"),
#         ("V1", "V3"),
#         ("V2", "V4"),
#         ("V4", "V5"),
#         ("V3", "V5"),
#         ("V4", "V6"),
#         ("V6", "V7"),
#         ("V5", "V7"),
#         ("V7", "V8"),
#         ("V6", "V8"),
#     ]
# )

# print(get_decycling_number(nt))
# print(get_decycling_number_mif(nt))
# print(get_decycling_number_mif_v2(nt))
# print(get_decycling_number_mif_v3(nt))

# graph = Network()
# graph.inherit_edge_colors(False)
# graph.options.edges.smooth.enabled = False
# graph.from_nx(nt)
# graph.get_node("V6")["color"] = "green"
#
# graph.show("graph.html", notebook=False)
