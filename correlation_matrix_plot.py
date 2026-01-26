#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt

def rates_from_confusion(TP, FP, FN, TN):
    FPR = FP / (FP + TN) if (FP + TN) > 0 else 0.0
    TPR = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    return FPR, TPR

def main():
    if len(sys.argv) < 5:
        print("Usage:")
        print("python3 plot_roc.py output.pdf TOOL1 TOOL2 ...")
        print("Then one line per tool: TP FP FN TN")
        sys.exit(1)

    out_pdf = sys.argv[1]
    name = sys.argv[2]

    # Split arguments into lines as pasted from shell
    args = sys.argv[3:]

    # First chunk: tool names (until we hit a number)
    tool_names = []
    i = 0
    while i < len(args) and not args[i][0].isdigit():
        tool_names.append(args[i])
        i += 1

    n_tools = len(tool_names)

    if n_tools == 0:
        print("Error: no tool names found")
        sys.exit(1)

    remaining = args[i:]

    if len(remaining) != 4 * n_tools:
        print("Error: number of TP FP FN TN values does not match number of tools")
        print(f"Expected {4*n_tools} numbers, got {len(remaining)}")
        sys.exit(1)

    plt.figure()

    idx = 0
    for tool in tool_names:
        TP = float(remaining[idx])
        FP = float(remaining[idx + 1])
        FN = float(remaining[idx + 2])
        TN = float(remaining[idx + 3])
        idx += 4

        fpr, tpr = rates_from_confusion(TP, FP, FN, TN)

        print(f"{tool}: FPR={fpr:.6g}, TPR={tpr:.6g}")
        plt.scatter(fpr, tpr, label=tool)

    # Random classifier diagonal
    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)

    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC {name}")
    plt.legend()
    plt.grid(True)

    plt.savefig(out_pdf)
    print(f"Saved ROC plot to {out_pdf}")

if __name__ == "__main__":
    main()
