import matplotlib.pyplot as plt
import numpy as np

# filter entries
#infeasibility = np.array([10, 20, 45, 70])
#unconstrained_merit = np.array([39, 36, 32, 21])
#infeasibility_upper_bound = 100
infeasibility = np.array([10, 29.51275, 133.1133])
unconstrained_merit = np.array([54.59781, 39.84375, 0])
infeasibility_upper_bound = 166.392

# filter parameters
infeasibility_margin = 0.95 # 0.999
envelope_slope = 0.01 # 0.001

# envelope calculation
margin_infeasibility = infeasibility_margin * infeasibility
outer_cusps = unconstrained_merit - envelope_slope*margin_infeasibility
margin_infeasibility = np.append(margin_infeasibility, infeasibility_upper_bound)
inner_cusps = unconstrained_merit - envelope_slope*margin_infeasibility[1:]

# colors
filter_borders_color = '#666666' # gray
filter_envelope_color = '#cccccc' # lightgray

# filter entries
plt.scatter(infeasibility, unconstrained_merit)
# infeasibility upper bound
plt.axvline(x=infeasibility_upper_bound, color='r', linestyle='-')
plt.gca().set_xlim(left=0.)

def merge_lists(lst1, lst2):
   return np.array([[i, j] for i, j in zip(lst1, lst2)]).ravel()

# filter borders and envelope
number_entries = len(infeasibility)
for i in range(number_entries):
   # filter borders
   if 0 < i:
      plt.vlines(infeasibility[i], unconstrained_merit[i], unconstrained_merit[i-1], colors=filter_borders_color)
   if i < number_entries-1:
      plt.hlines(unconstrained_merit[i], infeasibility[i], infeasibility[i+1], colors=filter_borders_color)
   
   # envelope (vertical)
   #if 0 < i:
   #   plt.vlines(infeasibility_margin*infeasibility[i], outer_cusps[i], inner_cusps[i-1], colors=filter_envelope_color)
   # envelope (slope): y <= unconstrained_merit[i] - envelope_slope * x
   #plt.plot([margin_infeasibility[i], margin_infeasibility[i+1]], [outer_cusps[i], inner_cusps[i]], color=filter_envelope_color)
   
plt.hlines(unconstrained_merit[-1], infeasibility[-1], infeasibility_upper_bound, colors=filter_borders_color)

# acceptable area
x = merge_lists(margin_infeasibility, margin_infeasibility)
x = np.insert(x, 0, 0.)
x = x[0:-1]
y1 = merge_lists(outer_cusps, inner_cusps)
(min_y, max_y) = plt.ylim()
y1 = np.insert(y1, 0, max_y)
y1 = np.insert(y1, 0, max_y)
plt.fill_between(x, y1=y1, y2=min_y, facecolor='g', alpha=0.3)

plt.xlabel("infeasibility")
plt.ylabel("unconstrained merit")
plt.show()
