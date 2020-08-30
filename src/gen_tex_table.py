import sys
import itertools

if len(sys.argv) != 2:
    print("Usage: python " + sys.argv[0] + " <filename>")
    exit(1)

best = 1e9
    
def convert_row(row):
    global best
    best = min(best, row[0])
    row[0] = str(round(1e-5*row[0], 2))
    if len(row) > 1:
        row[1] = str(round(row[1]))
        row[2] = str(round(row[2]))
    return row
    
content = []    
with open(sys.argv[1], "r") as in_file:
    content = [convert_row([float(i) for i in line.split(' ') if i.strip()]) for line in in_file]
    
#All experiments    
#parameters = ["1.", "2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała", "liniowa"]
#grids = ["$S_1$", "$_{0.5}$", "$N_{kon}$", "$N_{sr}$"]
#steps = ["0.5", "0.1", "0.02"]
#starts = ["$g_0$", "$g_3$", "$g_{0,0.4}$", "$g_{0,0.55}$", "$g_{3,0,0.5}$"]

#C1 test
#parameters = ["1."]
#backends = ["{\\it lm\\/}"]
#discrs = ["stała", "liniowa"]
#grids = ["$S_1$", "$S_{0.5}$", "$N_{kon}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_3$"]

#Discretization test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała", "liniowa"]
#grids = ["$S_{0.5}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_{0,0.4}$"]

# Grid test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała"]
#grids = ["$S_1$", "$S_{0.5}$", "$N_{sr}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_{0,0.4}$"]

# h test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała"]
#grids = ["$S_{0.5}$"]
#steps = ["0.5", "0.1", "0.02"]
#starts = ["$g_0$", "$g_{0,0.4}$"]

# start test
parameters = ["2."]
backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
discrs = ["stała"]
grids = ["$S_{0.5}$"]
steps = ["0.1"]
starts = ["$g_0$", "$g_3$", "$g_{0,0.4}$", "$g_{0,0.55}$", "$g_{3,0,0.5}$"]

res = '''\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}
\\hline
Parametry & algorytm & dyskret. & siatka & $h$ & start & wynik & \\#iteracji & \\#wywołań $\\hat{J}$ \\\\
\\hline
'''
cases = [case + (tuple(result)) for (case,result) in
         zip(itertools.product(parameters, backends, discrs, grids, steps, starts), content)]

for case in cases:
    res += " & ".join(case) + " \\\\\n\\hline\n"

res += '\\end{tabular}\n'
print(res)

print("Najlepszy: " + str(best));
