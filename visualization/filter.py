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
filter_borders_color = 'black' # gray
# filter_envelope_color = '#cccccc' # lightgray
entries_color = 'black'
upper_bound_color = 'r'
acceptable_area_color = 'g'

# filter entries
plt.scatter(infeasibility, unconstrained_merit, c=entries_color)
# infeasibility upper bound
plt.axvline(x=infeasibility_upper_bound, color=upper_bound_color, linestyle='-')
plt.gca().set_xlim(left=0.)

def merge_lists(lst1, lst2):
   return np.array([[i, j] for i, j in zip(lst1, lst2)]).ravel()

# filter borders
number_entries = len(infeasibility)
for i in range(number_entries):
   if 0 < i:
      plt.vlines(x=infeasibility[i], ymin=unconstrained_merit[i], ymax=unconstrained_merit[i-1], colors=filter_borders_color)
   if i < number_entries-1:
      plt.hlines(unconstrained_merit[i], infeasibility[i], infeasibility[i+1], colors=filter_borders_color)
plt.hlines(unconstrained_merit[-1], infeasibility[-1], infeasibility_upper_bound, colors=filter_borders_color)

(ymin, ymax) = plt.ylim()

# first vertical line
plt.vlines(x=infeasibility[0], ymin=unconstrained_merit[0], ymax=ymax, colors=filter_borders_color)

# acceptable area
x = merge_lists(margin_infeasibility, margin_infeasibility)
x = np.insert(x, 0, 0.)
x = x[0:-1]
y1 = merge_lists(outer_cusps, inner_cusps)
y1 = np.insert(y1, 0, ymax)
y1 = np.insert(y1, 0, ymax)
plt.fill_between(x, y1=y1, y2=ymin, facecolor=acceptable_area_color, alpha=0.3)

plt.xlabel("infeasibility")
plt.ylabel("unconstrained merit")
plt.show()
