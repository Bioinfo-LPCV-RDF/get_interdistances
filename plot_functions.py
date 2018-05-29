import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import numpy as np
import sys
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties


def divide(a, b):
    if b == 0:
        return np.nan
    else: 
        return a/b

def negative_sets_histo (Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command) :
	fig = plt.figure(1,figsize= (18,10))
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax1 = fig.add_subplot(1,4,1)
	ax1.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax1.axis([0, 3, 0, 2])
	ax1.bar(indexes1, rate , width , color = 'green')
	ax1.set_ylabel('$\sum DR$n_+, ER$n_+, IR$n_+ / $\sum DR$n_-, ER$n_-, IR$n_- , respectively', color = 'green', size = 16)
	for a,b in zip(relative_DR,relative_DR_neg) :
		ax1 = fig.add_subplot(1,4,2)
		ax1.set_xlabel("base pairs between direct repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('DR$n_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("DR$n_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("DR$n_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "D$Rn_+$ = DRn number in the bound set")
		plt.text(2, 4, "D$Rn_-$ = DRn number in the unbound set")
		plt.text(2, 3, "DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		plt.text(2, 6, sys.argv)
		l = plt.axhline(y = 1)
		ax1.legend()
		ax2.legend()
		
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		ax1 = fig.add_subplot(1,4,3)
		ax1.set_xlabel("base pairs between everted repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('E$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("E$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("E$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "ER$n_+$ = ERn number in the bound set")
		plt.text(2, 4, "ER$n_-$ = ERn number in the unbound set")
		plt.text(2, 3, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		ax1.legend()
		ax2.legend()
	for a,b in zip(relative_IR,relative_IR_neg) :	
		ax1 = fig.add_subplot(1,4,4)
		ax1.set_xlabel("base pairs between inverted repeats", size = 16)
		indexes1 = range(Interdistance_maxValue + 1)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
		ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
		ax2 = ax1.twinx()
		ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
		ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
		ax2.set_ylabel('I$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
		ybox1 = TextArea("I$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox3 = TextArea("I$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
		ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
		anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
				bbox_transform = ax2.transAxes, borderpad = 0.)
		ax2.add_artist(anchored_ybox)
		indexes1 = np.arange(Interdistance_maxValue + 1)
		plt.xticks(indexes1 + width * 0.5 , indexes1)
		plt.text(2, 5, "IR$n_+$ = DRn number in the bound set")
		plt.text(2, 4, "IR$n_-$ = DRn number in the unbound set")
		plt.text(2, 3, "IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1)
		ax1.legend()
		ax2.legend()
	
	plt.show()
		
def negative_sets_curve(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command) :
	fig = plt.figure(1,figsize= (18,10))
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax1 = fig.add_subplot(1,4,1)
	ax1.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 + width * 0.5, ('DR', 'ER', 'IR'))
	ax1.axis([0, 3, 0, 2])
	ax1.bar(indexes1, rate , width , color = 'green')
	ax1.set_ylabel('$\sum DR$n_+, ER$n_+, IR$n_+ / $\sum DR$n_-, ER$n_-, IR$n_- , respectively', color = 'green', size = 16)
	plt.text(2, 3, command)
	for a,b in zip(relative_DR,relative_DR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,2)
		ax1.set_xlabel("base pairs between direct repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2,label=r"")
		ax1.set_ylabel("DR$n_+ frequence / DR$n_- frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend()
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,3)
		ax1.set_xlabel("base pairs between everted repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2, label= r"threshold : "+str(c))
		ax1.set_ylabel("$ERn_+$ frequence / $ERn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
		plt.legend()
	for a,b in zip(relative_IR,relative_IR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		ax1 = fig.add_subplot(1,4,4)
		ax1.set_xlabel("base pairs between inverted repeats", size = 16)
		ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
		plt.plot(indexes1, map(divide, a, b), lw=2,label=r"")
		ax1.set_ylabel("$IRn_+$ frequence / $IRn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend()
	plt.legend()
	plt.show()
	
def negative_sets_points(Interdistance_maxValue,relative_DR,relative_DR_neg,relative_ER,relative_ER_neg,relative_IR,relative_IR_neg,threshold,rate,command,output,load=False) :	
        matplotlib.rcParams['font.sans-serif']=['Arial']
        matplotlib.rcParams['font.family']="sans-serif"
        SMALL_SIZE=5.4
        plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
	fig = plt.figure(1,figsize= (3.54,2.2))
	# fig.suptitle(command, fontsize = 11, fontweight='bold')
	gs = gridspec.GridSpec(1, 4,wspace = 0.25, width_ratios=[1, 2,2,2],left=0.1,right=0.97) 
	indexes1 = range(Interdistance_maxValue + 1)
	width = 1
	ax0 = plt.subplot(gs[0])
        ax0.spines["top"].set_linewidth(0.3)
        ax0.spines["bottom"].set_linewidth(0.3)
        ax0.spines["right"].set_linewidth(0.3)
        ax0.spines["left"].set_linewidth(0.3)
        ax0.tick_params(length=1.5,width=0.5,pad=1.5)
	ax0.set_xlabel("", size = 16)
	indexes1 = np.arange(3)
	plt.xticks(indexes1 , ('DR', 'ER', 'IR'))
	turn = 0
        if not load:
            points = ['lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black']
            ax0.set_color_cycle(['lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            points = ['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue']
            ax0.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
	maxi = []
	for a in rate :
		maxi.append(max(a))
	m = max(maxi)
	ax0.axis([-1, 3 , 0, m + 0.5])
	for a in rate :
		plt.plot(indexes1, a, points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
		turn = turn + 1
	ax0.set_ylabel("Absolute enrichment", fontsize = 5.4)
	turn = 0
	#plt.text(-3.6, -0.2, command)	
	maxi = []
	for a,b in zip(relative_DR,relative_DR_neg) :
		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_ER,relative_ER_neg) :
		maxi.append(max(map(divide, a, b)))
	for a,b in zip(relative_IR,relative_IR_neg) :
		maxi.append(max(map(divide, a, b)))
	m = max(maxi)		
	ax1 = plt.subplot(gs[1])
        ax1.spines["top"].set_linewidth(0.3)
        ax1.spines["bottom"].set_linewidth(0.3)
        ax1.spines["right"].set_linewidth(0.3)
        ax1.spines["left"].set_linewidth(0.3)
        ax1.tick_params(length=1.5,width=0.5,pad=1.5)
        ax1.yaxis.set_visible(False)
        if not load:
            ax1.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            ax1.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
	for a,b in zip(relative_DR,relative_DR_neg) :
		indexes1 = np.arange(Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(2,2,2)
		ax1.set_xlabel("bp between DRs", fontsize = 5.4)
		ax1.axis([-1, Interdistance_maxValue + 1, 0, m + 0.5])
		#plt.Circle((indexes1, map(divide, a, b)), color = 'g', radius=50, fill=True)
		plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2, markeredgewidth = 0.0)
		#plt.title('DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
		
		# t = plt.text(0.5,0.5, 'DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$',horizontalalignment='center', style='italic',verticalalignment='center', transform = ax1.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.text(2, 3.5, "DRn frequence = DRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1,linewidth=0.1)
		for i in range (0,Interdistance_maxValue + 1) :
                    if i%5 == 0 :
			plt.axvline(x = i, linewidth=0.3)
		turn = turn + 1
	turn = 0
	ax2 = plt.subplot(gs[2])
        ax2.spines["top"].set_linewidth(0.3)
        ax2.spines["bottom"].set_linewidth(0.3)
        ax2.spines["right"].set_linewidth(0.3)
        ax2.spines["left"].set_linewidth(0.3)
        ax2.tick_params(length=1.5,width=0.5,pad=1.5)        
        ax2.set_ylabel("Normalized enrichment", fontsize = 5.4, labelpad=60)
        if not load:
            ax2.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            ax2.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])

	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	    
		indexes1 = range(Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(1,7,4)
		ax2.set_xlabel("bp between ERs", fontsize = 5.4)
		ax2.axis([-1, Interdistance_maxValue + 1, 0, m + 0.5])
                try :
                    int(c)
                    plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, label = r"threshold : " + str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)
                except:
                    plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, label = str(c), alpha = 0.70, markersize = 2,markeredgewidth = 0.0)
		#plt.plot(indexes1, map(divide, a, b), points[turn], marker='o', alpha = 0.70, markersize = 13,markeredgewidth = 0.0)
		#plt.figlegend( d,('lightsalmon', 'tomato', 'r'), loc = 'lower center', ncol=5, labelspacing=0. )

		# ax2.set_ylabel("ERn+ frequence / ERn- frequence", size = 16)
		# t = plt.text(0.5,0.5, 'ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', horizontalalignment='center', style='italic',verticalalignment='center', transform = ax2.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.title('ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
		#plt.text(0, 3.5, "ERn frequence = ERn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$")
		l = plt.axhline(y = 1,linewidth=0.1)
		plt.legend(bbox_to_anchor=(0.5,0.5), ncol=(len(threshold)//2), mode = "expand", fontsize=5.4 ,borderaxespad = -15)
		turn = turn + 1
		for i in range (0,Interdistance_maxValue + 1) :
                    if i%5 == 0 :
			plt.axvline(x = i, linewidth=0.3)
	turn = 0
	ax3 = plt.subplot(gs[3])
        ax3.spines["top"].set_linewidth(0.3)
        ax3.spines["bottom"].set_linewidth(0.3)
        ax3.spines["right"].set_linewidth(0.3)
        ax3.spines["left"].set_linewidth(0.3)
        ax3.tick_params(length=1.5,width=0.5,pad=1.5)
        if not load:
            ax3.set_color_cycle([ 'lightsalmon', 'tomato', 'r', 'brown', 'maroon', 'black'])
        else:
            ax3.set_color_cycle(['lightsalmon', 'skyblue', 'r', 'royalblue', 'maroon', 'darkblue'])
	for a,b in zip(relative_IR,relative_IR_neg) :
		indexes1 = range(Interdistance_maxValue + 1)
		#ax1 = fig.add_subplot(1,7,6)
		ax3.set_xlabel("bp between IRs", fontsize = 5.4)
		ax3.axis([-1, Interdistance_maxValue + 1, 0, m + 0.5])
		plt.plot(indexes1, map(divide, a, b), points[turn], marker='o',linewidth=0.5, alpha = 0.70, markersize = 2,markeredgewidth = 0.0)

		# ax3.set_ylabel("IRn+ frequence / IRn- frequence", size = 16)
		# t = plt.text(0.5,0.5, 'IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', horizontalalignment='center', style='italic',verticalalignment='center', transform = ax3.transAxes)
		# t.set_bbox(dict(facecolor='white', alpha=0.1, edgecolor='black'))
		#plt.title('IRn frequence = IRn / $\sum_{i=0}^{i=20}$ DR$_i$ + ER$_i$ + IR$_i$', fontsize = 14)
		l = plt.axhline(y = 1,linewidth=0.1)
		turn = turn + 1
		for i in range (0,Interdistance_maxValue + 1):
                    if i%5 == 0 :
			plt.axvline(x = i, linewidth=0.3)
	plt.subplots_adjust(top=0.85)
#        plt.tight_layout(fig)
	fig.savefig(output)

