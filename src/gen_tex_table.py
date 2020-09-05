import sys
import itertools

if not (len(sys.argv) in [2,3]):
    print("Usage: python " + sys.argv[0] + " <filename> [<filename>]")
    exit(1)

#if len(sys.argv) > 2:
#    with open(sys.argv[2], "r") as in_file:
#        grads = [float(line) for line in in_file]
    
best = 1e9

#i = 0

def convert_row(row):
    global best
    global i
    best = min(best, row[0])
    row[0] = str(round(1e-5*row[0], 2))
    if len(row) > 1:
        row[1] = str(round(row[1]))
        row[2] = str(round(row[2]))
        if len(row) > 4:
            ratio = row[3] / row[4]
            if ratio < 1e-2:
                row[4] = "{:.1e}".format(ratio)
            else:
                row[4] = str(round(ratio, 3))
#        if len(sys.argv) > 2:
#            row.append(str(round(grads[i] / row[3], 1)))
#            i += 1
        row[3] = str(round(1e-5*row[3],2))
    return row
    
content = []    
with open(sys.argv[1], "r") as in_file:
    content = [convert_row([float(i) for i in line.split(' ') if i.strip()]) for line in in_file]
    
#All experiments    
#parameters = ["(CC)", "(DC)"]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["$P_0$", "$P_1$"]
#grids = ["$S_1$", "$_{0.5}$", "$N_{kon}$", "$N_{sr}$"]
#steps = ["0.5", "0.1", "0.02"]
#starts = ["$g_0$", "$g_3$", "$g_{0,0.4,42.5}$", "$g_{0,0.55,42.5}$", "$g_{0.07,0.59,48.2}$", "$\\mathfrak{g}$"]
#precisions = ["$10^{-6}$", "$10^{-9}$"]

#C1 test
#parameters = ["(CC)"]
#backends = ["{\\it lm\\/}"]
#discrs = ["$P_0$", "$P_1$"]
#grids = ["$S_1$", "$S_{0.5}$", "$N_{kon}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_3$"]
#precisions = ["$10^{-6}$", "$10^{-9}$"]

# start test
#parameters = ["(DC)"]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["$P_0$"]
#grids = ["$S_{0.5}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_3$", "$g_{0,0.55,42.5}$", "$g_{0.07,0.59,48.2}$"]
#precisions = ["$10^{-6}$", "$10^{-9}$"]

#Discretization test
#parameters = ["(DC)"]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["$P_1$"]
#grids = ["$S_{0.5}$"]
#steps = ["0.1"]
#starts = ["$g_{0.07,0.59,48.2}$"]
#precisions = ["$10^{-9}$"]

# Grid test
#parameters = ["(DC)"]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["$P_0$"]
#grids = ["$S_1$", "$N_{sr}$"]
#steps = ["0.1"]
#starts = ["$g_{0.07,0.59,48.2}$"]
#precisions = ["$10^{-9}$"]

# h test
parameters = ["(DC)"]
backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
discrs = ["$P_0$"]
grids = ["$S_{0.5}$"]
steps = ["0.5", "0.02"]
starts = ["$g_{0.07,0.59,48.2}$"]
precisions = ["$10^{-9}$"]

# precision test
parameters = ["(DC)"]
backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
discrs = ["$P_0$"]
grids = ["$S_{0.5}$"]
steps = ["0.1"]
starts = ["$g_{0.07,0.59,48.2}$"]
precisions = ["$10^{-6}$", "$10^{-9}$", "$10^{-11}$", "$10^{-13}$"]

res = '''\\begin{tabular}{|c|c|c|c|c|c|c||c|c|c|c|c|}
\\hline
Param. & algorytm & aproks. & siatka & $h$ & start & Tol & $\hat{J}$ & iter & $\\#\\hat{J}$ & $\\norm{G}_1$ & $\\frac{\\norm{G}_1}{\\norm{G_0}_1}$ \\\\
\\hline
'''
cases = [case + (tuple(result)) for (case,result) in
         zip(itertools.product(parameters, backends, discrs, grids, steps, starts, precisions), content)]

starts_order = {"$g_0$": 0, "$g_3$":1, "$g_{0,0.4,42.5}$":2, "$g_{0,0.55,42.5}$":3, "$g_{0.07,0.59,48.2}$":4, "$\\mathfrak{g}$":5}
cases.sort(key=lambda tup: starts_order[tup[5]])

for case in cases:
    res += " & ".join(case) + " \\\\\n\\hline\n"

res += '\\end{tabular}\n'
print(res)

print("Najlepszy: " + str(best));
