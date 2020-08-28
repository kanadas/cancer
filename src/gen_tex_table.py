import sys
import itertools

if len(sys.argv) != 1:
    print("Usage: python " + sys.argv[0] + " <filename>")

content = []    
with open(sys.argv[1], "r") as in_file:
    content = [float(i) for line in in_file for i in line.split(' ') if i.strip()]

parameters = ["1.", "2."]
backends = ["{\\it lm}", "{\\it sqp}"]
discrs = ["stała", "liniowa"]
grids = ["krok 1", "krok 0.5", "gęstsze końce", "gęstszy środek"]
steps = ["0.5", "0.1", "0.02"]
starts = ["0", "$g_max$", "$0.4\\cdot\\1_{[42.5,200]}$", "$0.55\\cdot\\1_{[42.5,200]}$"]
#
#def product(ar_list):
#    if not ar_list:
#        yield ()
#    else:
#        for a in ar_list[0]:
#            for prod in product(ar_list[1:]):
#                yield (a,)+prod
#
res = '''\\begin{tabular}{|c|c|c|c|c|c|c|}
\\hline
parametry & algorytm & dyskret. & siatka & h & start & wynik \\\\
\\hline
'''

best = {"1.": 100000000, "2.": 100000000}

cases = [case + (str(result),) for (case,result) in
         zip(itertools.product(parameters, backends, discrs, grids, steps, starts), content)]

for case in cases:
    best[case[0]] = min(best[case[0]], float(case[6]))
    res += " & ".join(case) + " \\\\\n\\hline\n"

res += '\\end{tabular}\n'
print(res)

print("\n\n")
print("Wyniki: " + str(best));
