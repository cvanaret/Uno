import matplotlib.pyplot as plt
import numpy as np

infeasibility = np.array([0.5, 20, 45, 70])
unconstrained_merit = np.array([39, 36, 32, 21])
infeasibility_upper_bound = 100
infeasibility_margin = 0.95 # 0.999
envelope_slope = 0.01 # 0.001

# envelope calculation
margin_infeasibility = infeasibility_margin * infeasibility
outer_cusps = unconstrained_merit - envelope_slope*margin_infeasibility
margin_infeasibility = np.append(margin_infeasibility, infeasibility_upper_bound)
inner_cusps = unconstrained_merit - envelope_slope*margin_infeasibility[1:]

print('outer cusps y: ', outer_cusps)
print('inner cusps y: ', inner_cusps)

# colors
filter_borders_color = '#666666' # gray
filter_envelope_color = '#cccccc' # lightgray

# filter entries
plt.scatter(infeasibility, unconstrained_merit)
# infeasibility upper bound
plt.axvline(x=infeasibility_upper_bound, color='r', linestyle='-')

# filter borders and envelope
number_entries = len(infeasibility)
for i in range(number_entries):
   # filter borders
   if 0 < i:
      # x = infeasibility[i]
      plt.vlines(infeasibility[i], unconstrained_merit[i], unconstrained_merit[i-1], colors=filter_borders_color)
   if i < number_entries-1:
      plt.hlines(unconstrained_merit[i], infeasibility[i], infeasibility[i+1], colors=filter_borders_color)
   
   # envelope
   if 0 < i:
      # x = infeasibility_margin*infeasibility[i]
      plt.vlines(infeasibility_margin*infeasibility[i], outer_cusps[i], inner_cusps[i-1], colors=filter_envelope_color)
      
   # envelope: y <= unconstrained_merit[i] - envelope_slope * x
   plt.plot([margin_infeasibility[i], margin_infeasibility[i+1]], [outer_cusps[i], inner_cusps[i]], color=filter_envelope_color)
   

plt.hlines(unconstrained_merit[-1], infeasibility[-1], infeasibility_upper_bound, colors=filter_borders_color)

plt.xlabel("infeasibility")
plt.ylabel("unconstrained merit")
plt.show()
