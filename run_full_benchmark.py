import datetime
import multiprocessing
import sys
import time
import os
from queue import Empty as QueueEmpty

import networkx as nx
import numpy as np

from naive import get_decycling_number_naive
from razgon import get_decycling_number_razgon
from fomin import get_decycling_number_fomin
from xiao import get_decycling_number_xiao
from bafna import approx_decycling_number_bafna
from bar_yehuda import approx_decycling_number_bar_yehuda
from stanojevic import approx_decycling_number_stanojevic


sys.setrecursionlimit(10000)


def parse_npz_file(file_path, has_dn=False):
    graphs = []
    try:
        with np.load(file_path, allow_pickle=True) as data:
            if "graphs" in data:
                matrices = data["graphs"].tolist()
                for i, matrix in enumerate(matrices):
                    G = nx.from_scipy_sparse_array(matrix)
                    if has_dn:
                        dn = data["decycling_numbers"][i]
                        G.graph["decycling_number"] = dn
                    graphs.append(G)
            else:
                print(f"Error: No 'graphs' key found in '{file_path}'.")
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    return graphs


def run_method_get_time_and_result(method_func, graph_data, result_queue):
    try:
        start_time = time.perf_counter()
        result_value = method_func(graph_data)
        end_time = time.perf_counter()
        result_queue.put((end_time - start_time, result_value))
    except Exception as e:
        result_queue.put(e)


def run_method_get_result(method_func, graph_data, result_queue):
    try:
        result = method_func(graph_data)
        result_queue.put(result)
    except Exception as e:
        result_queue.put(e)


def run_approx_method_get_result_and_time(method_func, graph_data, result_queue):
    try:
        start_time = time.perf_counter()
        k = method_func(graph_data)
        end_time = time.perf_counter()
        duration = end_time - start_time
        result_queue.put((k, duration))
    except Exception as e:
        result_queue.put(e)


def benchmark_exact_exec_time(
    directory_path,
    methods_list,
    timeout_seconds,
    output_filename,
):

    benchmark_start_time = time.perf_counter()

    if not os.path.isdir(directory_path):
        print(f"Error : Folder '{directory_path}' does not exist.")
        return

    files_to_process = sorted(
        f for f in os.listdir(directory_path) if f.endswith(".npz")
    )

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
                file_results_values = []
                file_header_line = f"Treating: {filename}"
                file_separator = "=" * 80

                f_out.write(f"{file_separator}\n")
                f_out.write(f"{file_header_line}\n")
                f_out.write(f"{file_separator}\n\n")

                filepath = os.path.join(directory_path, filename)

                try:
                    graphs_in_file = parse_npz_file(filepath)
                    total_graphs = len(graphs_in_file)
                    if total_graphs == 0:
                        # Skip le fichier
                        continue
                    f_out.write(f"  File loaded: {total_graphs} graphs.\n")

                except Exception as e:
                    print(f"Error loading graphs from '{filename}': {e}")
                    continue

                results = {
                    m.__name__: {"times": [], "timeouts": 0, "errors": 0, "success": 0}
                    for m in methods_list
                }

                for i, graph in enumerate(graphs_in_file):
                    current_graph_result_value = None
                    for method in methods_list:
                        method_name = method.__name__
                        run_times = []
                        timeout_occurred = False
                        error_occurred = False

                        for _ in range(3):
                            result_queue = multiprocessing.Queue()
                            p = multiprocessing.Process(
                                target=run_method_get_time_and_result,
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
                                        duration, val = result
                                        run_times.append(duration)
                                        if current_graph_result_value is None:
                                            current_graph_result_value = val
                                        elif current_graph_result_value != val:
                                            print(
                                                f"\n Logic error: Inconsistent results for graph {i + 1} in file '{filename}' between methods. "
                                            )
                                            print(
                                                f"Method {method_name} returned {val}, previous was {current_graph_result_value}."
                                            )
                                            return
                                except QueueEmpty:
                                    error_occurred = True
                                    break

                        if run_times:
                            results[method_name]["times"].extend(run_times)
                        elif timeout_occurred:
                            results[method_name]["timeouts"] += 1
                        elif error_occurred:
                            results[method_name]["errors"] += 1

                    if current_graph_result_value is None:
                        file_results_values.append("Error/Timeout")
                    else:
                        file_results_values.append(current_graph_result_value)

                f_out.write("\n  Finished tests. See results.\n\n")

                col_methode = max(len(m.__name__) for m in methods_list) + 2
                col_methode = max(col_methode, 10)
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

                f_out.write(table_header + "\n")
                f_out.write(table_separator + "\n")

                for method_name, data in results.items():
                    times_sec = np.array(data["times"])
                    num_success = len(times_sec)
                    num_timeouts = data["timeouts"]
                    num_errors = data["errors"]

                    timeout_percent = (num_timeouts / total_graphs) * 100
                    error_percent = (num_errors / total_graphs) * 100
                    success_percent = ((num_success / 3) / total_graphs) * 100

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

                    f_out.write(row_data + "\n")
                f_out.write("\n")

            benchmark_end_time = time.perf_counter()
            total_seconds = benchmark_end_time - benchmark_start_time
            formatted_time = str(datetime.timedelta(seconds=total_seconds))
            f_out.write(f"\n{file_separator}\n")
            f_out.write("Benchmark over.\n")
            f_out.write(f"Total execution time: {formatted_time}\n")

    except IOError as e:
        print(f"Error: Impossible to write in '{output_filename}': {e}")
    except Exception as e:
        print(f"Error : {e}")


def benchmark_approximation_quality(
    directory_path,
    approx_methods_list,
    exact_method_func,
    timeout_seconds,
    output_filename,
):
    benchmark_start_time = time.perf_counter()

    if not os.path.isdir(directory_path):
        print(f"Error : Folder '{directory_path}' does not exist.")
        return

    files_to_process = [f for f in os.listdir(directory_path) if f.endswith(".npz")]

    if not files_to_process:
        print(f"No file to process has been found in folder '{directory_path}'.")
        return

    try:
        with open(output_filename, "w", encoding="utf-8") as f_out:
            f_out.write("=" * 100 + "\n")
            f_out.write(f"APPROXIMATION QUALITY BENCHMARK\n")
            f_out.write(f"Date: {time.ctime()}\n")
            f_out.write(f"Ref Exact Method: {exact_method_func.__name__}\n")
            f_out.write(f"Timeout per graph: {timeout_seconds}s\n")
            f_out.write("=" * 100 + "\n\n")

            for filename in files_to_process:
                file_header_line = f"Treating: {filename}"
                file_separator = "=" * 100
                f_out.write(f"{file_separator}\n")
                f_out.write(f"{file_header_line}\n")
                f_out.write(f"{file_separator}\n\n")
                filepath = os.path.join(directory_path, filename)

                try:
                    graphs_in_file = parse_npz_file(filepath)
                    total_graphs = len(graphs_in_file)
                    if total_graphs == 0:
                        continue
                    f_out.write(f"  File loaded: {total_graphs} graphs.\n")
                except Exception as e:
                    print(f"Error loading graphs from '{filename}': {e}")
                    continue

                ratios_results = {m.__name__: [] for m in approx_methods_list}
                times_results = {m.__name__: [] for m in approx_methods_list}
                meta_stats = {
                    m.__name__: {"timeouts": 0, "errors": 0, "success": 0}
                    for m in approx_methods_list
                }
                skipped_graphs_exact_fail = 0

                for i, graph in enumerate(graphs_in_file):
                    exact_queue = multiprocessing.Queue()
                    p_exact = multiprocessing.Process(
                        target=run_method_get_result,
                        args=(exact_method_func, graph, exact_queue),
                    )
                    p_exact.start()
                    p_exact.join(timeout=timeout_seconds)

                    k_opt = None
                    if p_exact.is_alive():
                        p_exact.terminate()
                        p_exact.join()
                        skipped_graphs_exact_fail += 1
                        # Skip ce graphe
                        continue

                    try:
                        res_exact = exact_queue.get_nowait()
                        if isinstance(res_exact, Exception):
                            skipped_graphs_exact_fail += 1
                            continue
                        k_opt = res_exact
                    except QueueEmpty:
                        skipped_graphs_exact_fail += 1
                        continue

                    for method in approx_methods_list:
                        method_name = method.__name__
                        best_k_approx = float("inf")
                        run_times_per_graph = []
                        valid_run_found = False
                        timeout_occurred = False
                        error_occurred = False

                        for _ in range(3):
                            approx_queue = multiprocessing.Queue()
                            p_approx = multiprocessing.Process(
                                target=run_approx_method_get_result_and_time,
                                args=(method, graph, approx_queue),
                            )
                            p_approx.start()
                            p_approx.join(timeout=timeout_seconds)

                            if p_approx.is_alive():
                                p_approx.terminate()
                                p_approx.join()
                                timeout_occurred = True
                                break
                            elif p_approx.exitcode != 0:
                                error_occurred = True
                                break
                            else:
                                try:
                                    res = approx_queue.get_nowait()
                                    if isinstance(res, Exception):
                                        error_occurred = True
                                        break
                                    else:
                                        res_k, res_t = res
                                        valid_run_found = True
                                        run_times_per_graph.append(res_t)
                                        if res_k < best_k_approx:
                                            best_k_approx = res_k
                                except QueueEmpty:
                                    error_occurred = True
                                    break

                        if valid_run_found:
                            meta_stats[method_name]["success"] += 1
                            if k_opt == 0:
                                ratio = 1.0 if best_k_approx == 0 else 999.0
                            else:
                                ratio = best_k_approx / k_opt
                            ratios_results[method_name].append(ratio)
                            if run_times_per_graph:
                                times_results[method_name].extend(run_times_per_graph)

                        else:
                            if timeout_occurred:
                                meta_stats[method_name]["timeouts"] += 1
                            else:
                                meta_stats[method_name]["errors"] += 1

                if skipped_graphs_exact_fail > 0:
                    print(
                        f"  WARNING: {skipped_graphs_exact_fail} graphs skipped because Exact Method failed/timeout."
                    )
                    f_out.write(
                        f"  WARNING: {skipped_graphs_exact_fail} graphs skipped because Exact Method failed/timeout.\n"
                    )

                f_out.write("\n Finished tests. See results.\n\n")

                col_methode = max(len(m.__name__) for m in approx_methods_list) + 2
                col_methode = max(col_methode, 15)
                col_val = 14
                col_pct = 12
                col_exact = 10
                col_time = 14

                table_header = (
                    f"| {'Method':<{col_methode}} "
                    f"| {'Exact Matches':>{col_exact}} "
                    f"| {'Min Ratio':>{col_val}} "
                    f"| {'Max Ratio':>{col_val}} "
                    f"| {'Mean Ratio':>{col_val}} "
                    f"| {'Std Ratio':>{col_val}} "
                    f"| {'Median Ratio':>{col_val}} "
                    f"| {'Mean Time (s)':>{col_time}} "
                    f"| {'Median Time (s)':>{col_time}} "
                    f"| {'Success (%)':>{col_pct}} "
                    f"| {'Timeout (%)':>{col_pct}} |"
                )

                table_separator = (
                    f"|{'-' * (col_methode + 1)}"
                    f"|{'-' * (col_exact + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_pct + 1)}|"
                )

                f_out.write(table_header + "\n")
                f_out.write(table_separator + "\n")

                effective_total = total_graphs - skipped_graphs_exact_fail
                if effective_total <= 0:
                    effective_total = 1

                for method_name, ratios_list in ratios_results.items():
                    ratios_arr = np.array(ratios_list)
                    times_list = times_results.get(method_name, [])
                    times_arr = np.array(times_list)
                    stats = meta_stats[method_name]
                    success_pct = (stats["success"] / effective_total) * 100
                    timeout_pct = (stats["timeouts"] / effective_total) * 100

                    if len(ratios_arr) > 0:
                        r_min = f"{np.min(ratios_arr):.4f}"
                        r_max = f"{np.max(ratios_arr):.4f}"
                        r_mean = f"{np.mean(ratios_arr):.4f}"
                        r_median = f"{np.median(ratios_arr):.4f}"
                        r_std = f"{np.std(ratios_arr):.4f}"
                        r_exact_matches = np.sum(ratios_arr == 1.0)
                    else:
                        r_mean = r_median = r_max = r_min = r_std = r_exact_matches = (
                            "---"
                        )

                    if len(times_arr) > 0:
                        rt_mean = f"{np.mean(times_arr):.6f}"
                        rt_median = f"{np.median(times_arr):.6f}"
                    else:
                        rt_mean = rt_median = "---"

                    row_data = (
                        f"| {method_name:<{col_methode}} "
                        f"| {r_exact_matches:>{col_exact}} "
                        f"| {r_min:>{col_val}} "
                        f"| {r_max:>{col_val}} "
                        f"| {r_mean:>{col_val}} "
                        f"| {r_std:>{col_val}} "
                        f"| {r_median:>{col_val}} "
                        f"| {rt_mean:>{col_time}} "
                        f"| {rt_median:>{col_time}} "
                        f"| {success_pct:>{col_pct - 3}.2f} % "
                        f"| {timeout_pct:>{col_pct - 3}.2f} % |"
                    )

                    f_out.write(row_data + "\n")

                f_out.write("\n")

            benchmark_end_time = time.perf_counter()
            total_seconds = benchmark_end_time - benchmark_start_time
            formatted_time = str(datetime.timedelta(seconds=total_seconds))
            end_msg = (
                f"\n{file_separator}\n"
                f"Quality Benchmark over.\n"
                f"Total execution time: {formatted_time}\n"
                f"{file_separator}"
            )
            f_out.write(end_msg + "\n")

    except IOError as e:
        print(f"Error: Impossible to write in '{output_filename}': {e}")
    except Exception as e:
        print(f"Error global : {e}")


def benchmark_approximation_quality_with_dn(
    directory_path,
    approx_methods_list,
    timeout_seconds,
    output_filename,
):
    benchmark_start_time = time.perf_counter()

    if not os.path.isdir(directory_path):
        print(f"Error : Folder '{directory_path}' does not exist.")
        return

    files_to_process = [f for f in os.listdir(directory_path) if f.endswith(".npz")]
    if not files_to_process:
        print(f"No file to process has been found in folder '{directory_path}'.")
        return

    try:
        with open(output_filename, "w", encoding="utf-8") as f_out:
            f_out.write("=" * 100 + "\n")
            f_out.write(f"APPROXIMATION QUALITY BENCHMARK\n")
            f_out.write(f"Date: {time.ctime()}\n")
            f_out.write(f"Timeout per graph: {timeout_seconds}s\n")
            f_out.write("=" * 100 + "\n\n")

            for filename in files_to_process:
                file_header_line = f"Treating: {filename}"
                file_separator = "=" * 100
                f_out.write(f"{file_separator}\n")
                f_out.write(f"{file_header_line}\n")
                f_out.write(f"{file_separator}\n\n")
                filepath = os.path.join(directory_path, filename)

                try:
                    graphs_in_file = parse_npz_file(filepath, has_dn=True)
                    total_graphs = len(graphs_in_file)
                    if total_graphs == 0:
                        continue
                    f_out.write(f"  File loaded: {total_graphs} graphs.\n")
                except Exception as e:
                    print(f"Error loading graphs from '{filename}': {e}")
                    continue

                ratios_results = {m.__name__: [] for m in approx_methods_list}
                times_results = {m.__name__: [] for m in approx_methods_list}
                meta_stats = {
                    m.__name__: {"timeouts": 0, "errors": 0, "success": 0}
                    for m in approx_methods_list
                }
                skipped_graphs_exact_fail = 0

                for i, graph in enumerate(graphs_in_file):
                    k_opt = graph.graph.get("decycling_number")

                    for method in approx_methods_list:
                        method_name = method.__name__
                        best_k_approx = float("inf")
                        run_times_per_graph = []
                        valid_run_found = False
                        timeout_occurred = False
                        error_occurred = False

                        for _ in range(3):
                            approx_queue = multiprocessing.Queue()
                            p_approx = multiprocessing.Process(
                                target=run_approx_method_get_result_and_time,
                                args=(method, graph, approx_queue),
                            )
                            p_approx.start()
                            p_approx.join(timeout=timeout_seconds)

                            if p_approx.is_alive():
                                p_approx.terminate()
                                p_approx.join()
                                timeout_occurred = True
                                break
                            elif p_approx.exitcode != 0:
                                error_occurred = True
                                break
                            else:
                                try:
                                    res = approx_queue.get_nowait()
                                    if isinstance(res, Exception):
                                        error_occurred = True
                                        break
                                    else:
                                        res_k, res_t = res
                                        valid_run_found = True
                                        run_times_per_graph.append(res_t)
                                        if res_k < best_k_approx:
                                            best_k_approx = res_k
                                except QueueEmpty:
                                    error_occurred = True
                                    break

                        if valid_run_found:
                            meta_stats[method_name]["success"] += 1
                            if k_opt == 0:
                                ratio = 1.0 if best_k_approx == 0 else 999.0
                            else:
                                ratio = best_k_approx / k_opt

                            ratios_results[method_name].append(ratio)
                            if run_times_per_graph:
                                times_results[method_name].extend(run_times_per_graph)

                        else:
                            if timeout_occurred:
                                meta_stats[method_name]["timeouts"] += 1
                            else:
                                meta_stats[method_name]["errors"] += 1

                if skipped_graphs_exact_fail > 0:
                    print(
                        f"  WARNING: {skipped_graphs_exact_fail} graphs skipped because Exact Method failed/timeout."
                    )
                    f_out.write(
                        f"  WARNING: {skipped_graphs_exact_fail} graphs skipped because Exact Method failed/timeout.\n"
                    )

                f_out.write("\n Finished tests. See results.\n\n")

                col_methode = max(len(m.__name__) for m in approx_methods_list) + 2
                col_methode = max(col_methode, 15)
                col_val = 14
                col_pct = 12
                col_exact = 10
                col_time = 14

                table_header = (
                    f"| {'Method':<{col_methode}} "
                    f"| {'Exact Matches':>{col_exact}} "
                    f"| {'Min Ratio':>{col_val}} "
                    f"| {'Max Ratio':>{col_val}} "
                    f"| {'Mean Ratio':>{col_val}} "
                    f"| {'Std Ratio':>{col_val}} "
                    f"| {'Median Ratio':>{col_val}} "
                    f"| {'Mean Time (s)':>{col_time}} "
                    f"| {'Median Time (s)':>{col_time}} "
                    f"| {'Success (%)':>{col_pct}} "
                    f"| {'Timeout (%)':>{col_pct}} |"
                )

                table_separator = (
                    f"|{'-' * (col_methode + 1)}"
                    f"|{'-' * (col_exact + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_val + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_time + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_pct + 1)}|"
                )

                f_out.write(table_header + "\n")
                f_out.write(table_separator + "\n")

                effective_total = total_graphs - skipped_graphs_exact_fail
                if effective_total <= 0:
                    effective_total = 1

                for method_name, ratios_list in ratios_results.items():
                    ratios_arr = np.array(ratios_list)
                    times_list = times_results.get(method_name, [])
                    times_arr = np.array(times_list)
                    stats = meta_stats[method_name]
                    success_pct = (stats["success"] / effective_total) * 100
                    timeout_pct = (stats["timeouts"] / effective_total) * 100

                    if len(ratios_arr) > 0:
                        r_min = f"{np.min(ratios_arr):.4f}"
                        r_max = f"{np.max(ratios_arr):.4f}"
                        r_mean = f"{np.mean(ratios_arr):.4f}"
                        r_median = f"{np.median(ratios_arr):.4f}"
                        r_std = f"{np.std(ratios_arr):.4f}"
                        r_exact_matches = np.sum(ratios_arr == 1.0)
                    else:
                        r_mean = r_median = r_max = r_min = r_std = r_exact_matches = (
                            "---"
                        )

                    if len(times_arr) > 0:
                        rt_mean = f"{np.mean(times_arr):.6f}"
                        rt_median = f"{np.median(times_arr):.6f}"
                    else:
                        rt_mean = rt_median = "---"

                    row_data = (
                        f"| {method_name:<{col_methode}} "
                        f"| {r_exact_matches:>{col_exact}} "
                        f"| {r_min:>{col_val}} "
                        f"| {r_max:>{col_val}} "
                        f"| {r_mean:>{col_val}} "
                        f"| {r_std:>{col_val}} "
                        f"| {r_median:>{col_val}} "
                        f"| {rt_mean:>{col_time}} "
                        f"| {rt_median:>{col_time}} "
                        f"| {success_pct:>{col_pct - 3}.2f} % "
                        f"| {timeout_pct:>{col_pct - 3}.2f} % |"
                    )

                    f_out.write(row_data + "\n")

                f_out.write("\n")

            benchmark_end_time = time.perf_counter()
            total_seconds = benchmark_end_time - benchmark_start_time
            formatted_time = str(datetime.timedelta(seconds=total_seconds))
            end_msg = (
                f"\n{file_separator}\n"
                f"Quality Benchmark over.\n"
                f"Total execution time: {formatted_time}\n"
                f"{file_separator}"
            )
            f_out.write(end_msg + "\n")

    except IOError as e:
        print(f"Error: Impossible to write in '{output_filename}': {e}")
    except Exception as e:
        print(f"Error global : {e}")


def benchmark_approx_comparison(
    directory_path,
    methods_list,
    timeout_seconds,
    output_filename="benchmark_comparison_results.txt",
):

    benchmark_start_time = time.perf_counter()

    if not os.path.isdir(directory_path):
        print(f"Error : Folder '{directory_path}' does not exist.")
        return

    files_to_process = [f for f in os.listdir(directory_path) if f.endswith(".npz")]
    if not files_to_process:
        print(f"No file to process has been found in folder '{directory_path}'.")
        return

    try:
        with open(output_filename, "w", encoding="utf-8") as f_out:
            f_out.write("=" * 100 + "\n")
            f_out.write(f"APPROXIMATION COMPARISON & RANKING\n")
            f_out.write(f"Date: {time.ctime()}\n")
            f_out.write(
                f"Ranking logic: Based on frequency of finding the MINIMUM result among competitors.\n"
            )
            f_out.write(f"Timeout per graph: {timeout_seconds}s\n")
            f_out.write("=" * 100 + "\n\n")

            for filename in files_to_process:
                file_header_line = f"Treating: {filename}"
                file_separator = "=" * 100
                f_out.write(f"{file_separator}\n")
                f_out.write(f"{file_header_line}\n")
                f_out.write(f"{file_separator}\n\n")
                filepath = os.path.join(directory_path, filename)

                try:
                    graphs_in_file = parse_npz_file(filepath)
                    total_graphs = len(graphs_in_file)
                    if total_graphs == 0:
                        continue
                    f_out.write(f"  File loaded: {total_graphs} graphs.\n")
                except Exception as e:
                    print(f"Error loading graphs from '{filename}': {e}")
                    continue

                stats = {
                    m.__name__: {
                        "times": [],
                        "results": [],
                        "wins": 0,
                        "timeouts": 0,
                        "errors": 0,
                        "success": 0,
                    }
                    for m in methods_list
                }

                for i, graph in enumerate(graphs_in_file):
                    current_graph_results = {}

                    for method in methods_list:
                        method_name = method.__name__

                        best_run_k = float("inf")
                        run_times_per_graph = []

                        run_valid = False
                        timeout_occurred = False
                        error_occurred = False

                        for _ in range(3):
                            q = multiprocessing.Queue()
                            p = multiprocessing.Process(
                                target=run_approx_method_get_result_and_time,
                                args=(method, graph, q),
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
                                    res = q.get_nowait()
                                    if isinstance(res, Exception):
                                        error_occurred = True
                                        break
                                    else:
                                        k_val, t_val = res
                                        run_times_per_graph.append(t_val)
                                        if k_val < best_run_k:
                                            best_run_k = k_val
                                            run_valid = True
                                        elif k_val == best_run_k:
                                            run_valid = True
                                except QueueEmpty:
                                    error_occurred = True
                                    break

                        if run_valid:
                            stats[method_name]["times"].extend(run_times_per_graph)
                            stats[method_name]["results"].append(best_run_k)
                            stats[method_name]["success"] += 1
                            current_graph_results[method_name] = best_run_k
                        else:
                            if timeout_occurred:
                                stats[method_name]["timeouts"] += 1
                            else:
                                stats[method_name]["errors"] += 1

                    if current_graph_results:
                        min_k_global = min(current_graph_results.values())
                        for m_name, k_val in current_graph_results.items():
                            if k_val == min_k_global:
                                stats[m_name]["wins"] += 1

                f_out.write("\n  Results table (Sorted by Best Solution Rate):\n\n")

                col_methode = max(len(m.__name__) for m in methods_list) + 2
                col_methode = max(col_methode, 15)
                col_time = 12
                col_pct = 14

                table_header = (
                    f"| {'Method':<{col_methode}} "
                    f"| {'Best Sol %':>{col_pct}} "
                    f"| {'Mean Time (s)':>{col_pct}} "
                    f"| {'Median Time (s)':>{col_pct}} "
                    f"| {'Success (%)':>{col_time}} |"
                )
                table_separator = (
                    f"|{'-' * (col_methode + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_pct + 1)}"
                    f"|{'-' * (col_time + 1)}|"
                )

                rows = []
                for method_name, data in stats.items():
                    num_runs = len(data["times"])
                    num_success = data["success"]
                    success_pct = (num_success / total_graphs) * 100
                    win_pct = (data["wins"] / total_graphs) * 100
                    times_arr = np.array(data["times"])

                    if num_runs > 0:
                        s_mean = np.mean(times_arr)
                        s_median = np.median(times_arr)
                    else:
                        s_mean = float("inf")
                        s_median = float("inf")

                    row_obj = {
                        "name": method_name,
                        "win_pct": win_pct,
                        "mean_time": s_mean,
                        "median_time": s_median,
                        "success_pct": success_pct,
                    }
                    rows.append(row_obj)

                rows.sort(key=lambda x: (-x["win_pct"], x["mean_time"]))
                f_out.write(table_header + "\n")
                f_out.write(table_separator + "\n")

                for row in rows:
                    mean_str = (
                        f"{row['mean_time']:.6f}"
                        if row["mean_time"] != float("inf")
                        else "---"
                    )
                    med_str = (
                        f"{row['median_time']:.6f}"
                        if row["median_time"] != float("inf")
                        else "---"
                    )
                    line = (
                        f"| {row['name']:<{col_methode}} "
                        f"| {row['win_pct']:>{col_pct - 3}.2f} % "
                        f"| {mean_str:>{col_pct}} "
                        f"| {med_str:>{col_pct}} "
                        f"| {row['success_pct']:>{col_time - 3}.2f} % |"
                    )
                    f_out.write(line + "\n")

                f_out.write("\n")

            benchmark_end_time = time.perf_counter()
            total_seconds = benchmark_end_time - benchmark_start_time
            formatted_time = str(datetime.timedelta(seconds=total_seconds))
            end_msg = (
                f"\n{file_separator}\n"
                f"Comparison Benchmark over.\n"
                f"Total execution time: {formatted_time}\n"
                f"{file_separator}"
            )
            f_out.write(end_msg + "\n")

    except IOError as e:
        print(f"Error: Impossible to write in '{output_filename}': {e}")
    except Exception as e:
        print(f"Error global : {e}")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    TIMEOUT_MAX = 600

    EXACT_METHODS = [
        get_decycling_number_razgon,
        get_decycling_number_fomin,
        get_decycling_number_xiao,
    ]

    APPROX_METHODS = [
        approx_decycling_number_bafna,
        approx_decycling_number_bar_yehuda,
        approx_decycling_number_stanojevic,
    ]

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/small for naive",
        methods_list=[get_decycling_number_naive] + EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_small_for_naive.txt",
    )

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/random graphs density",
        methods_list=EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_random_density.txt",
    )

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/density",
        methods_list=EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_density.txt",
    )

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/diameter",
        methods_list=EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_diameter.txt",
    )

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/radius",
        methods_list=EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_radius.txt",
    )

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/treewidth",
        methods_list=EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_treewidth.txt",
    )

    benchmark_exact_exec_time(
        directory_path="Benchmark graphs/vertex connectivity",
        methods_list=EXACT_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_vert_conn.txt",
    )

    benchmark_approximation_quality(
        directory_path="Benchmark graphs/random graphs density",
        approx_methods_list=APPROX_METHODS,
        exact_method_func=get_decycling_number_xiao,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_approx_random_density.txt",
    )

    benchmark_approximation_quality(
        directory_path="Benchmark graphs/density",
        approx_methods_list=APPROX_METHODS,
        exact_method_func=get_decycling_number_xiao,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_approx_density.txt",
    )

    benchmark_approximation_quality_with_dn(
        directory_path="Benchmark graphs/more vertices with dn",
        approx_methods_list=APPROX_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_approx_more_v_dn.txt",
    )

    benchmark_approx_comparison(
        directory_path="Benchmark graphs/random big",
        methods_list=APPROX_METHODS,
        timeout_seconds=TIMEOUT_MAX,
        output_filename="Benchmark results/final_ben_approx_random_big.txt",
    )
