#CC test
#parameters = ["1."]
#backends = ["{\\it lm\\/}"]
#discrs = ["stała", "liniowa"]
#grids = ["$S_1$", "$S_{0.5}$", "$N_{kon}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_3$"]
cc_del = 2*3*2

#Discretization test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała", "liniowa"]
#grids = ["$S_{0.5}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_{0,0.4}$"]
discr_del = 2*2

# Grid test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała"]
#grids = ["$S_1$", "$S_{0.5}$", "$N_{sr}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_{0,0.4}$"]
grid_del = 3*2

# h test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała"]
#grids = ["$S_{0.5}$"]
#steps = ["0.5", "0.1", "0.02"]
#starts = ["$g_0$", "$g_{0,0.4}$"]
h_del = 3*2

# start test
#parameters = ["2."]
#backends = ["{\\it lm\\/}", "{\\it sqp\\/}"]
#discrs = ["stała"]
#grids = ["$S_{0.5}$"]
#steps = ["0.1"]
#starts = ["$g_0$", "$g_3$", "$g_{0,0.4}$", "$g_{0,0.55}$", "$g_{3,0,0.5}$"]
start_del = 5

with open("res/res_start_lm", "r") as in_file:
    new_data = in_file.readlines()

with open("res/result_start", "r") as out_file:
    old_data = out_file.readlines()[start_del:]
    
with open("res/result_start", "w") as out_file:
    out_file.writelines(new_data + old_data)

