from matplotlib import pyplot as plt

plotl5_5 = [0.070031,0.105487,0.083339]
plotl5_9 = [0.450023,0.066627,0.463515]
plotl2_5 = [0.989226,0.718715,4.594056]
plotl2_9 = [4.794911,1.001329,5.452591]
#plotl5_5 = eval(input())
#plotl5_9 = eval(input())
#plotl2_5 = eval(input())
#plotl2_9 = eval(input())

plt.xlabel("Configuration")
plt.ylabel("Time in microseconds")

labels=["512^2 elements\nwith stencil 5","512^2 elements\nwith stencil 9","2048^2 elements\nwith stencil 5","2048^2 elements\nwith stencil 9"]
plt.boxplot([plotl5_5, plotl5_9, plotl2_5, plotl2_9],labels=labels)
plt.legend()
plt.yscale("log")
plt.title("[Stencil,Data points] vs Time")
plt.show()
plt.savefig("savedplot.png", format="png")
