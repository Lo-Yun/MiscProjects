import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure()
fig_row = 4
fig_col = 2
gs = gridspec.GridSpec(fig_row, fig_col)  ## 3 rows & 2 columns

set_energy_axis = [-10, 170, 3150, 3380]  ## format: [x_min, x_max, y_min, y_max]
set_state_ax = [-10, 170, 0, 1]    ## format: [x_min, x_max, y_min, y_max]
max_x = 160
xinterval = 10**(-2)

right_ylabel = True

################
## Sizes control
title_size = 8
ylabel_size = 8
##tick_size = 6
xtick_size = 7
y1tick_size = 7
y2tick_size = 8

for x_cell in (0, 1):
    ## Initialize
    x_plot = []
    xf_plot, xm_plot, xb_plot = [], [], []  # front-end, middle, back-end

    ev0 = []
    ev0f, ev0m, ev0b = [], [], []  # front-end, middle, back-end

    ev1 = []
    ev1f, ev1m, ev1b = [], [], []  # front-end, middle, back-end

    ev2 = []
    ev2f, ev2m, ev2b = [], [], []  # front-end, middle, back-end

    straight_Ea, straight_Eb, straight_Ec = [], [], []

    high_state_b0, mid_state_b0, low_state_b0 = [], [], []
    high_state_b1, mid_state_b1, low_state_b1 = [], [], []
    high_state_b2, mid_state_b2, low_state_b2 = [], [], []
    #state = {'high':{'wght1':[], }}
    
    if x_cell == 0:  ## at the begining of 2x2 region
        ####================####
        #### The 2x2 Matrix ####
        H = np.zeros([2,2])
        E0a = 3340
        E0b = 3220
        H[0,1], H[1,0] = 35, 35

        #############################################
        ## Update the matrix elements, then solve it!
        for x in np.arange( 0, max_x, xinterval ):
            H[0][0] = E0a - x
            H[1][1] = E0b + 0.5 * x
            straight_Ea.append(H[0][0])
            straight_Eb.append(H[1][1])
            x_plot.append(x)
            ev, ef = linalg.eig(H)
            #print 'x = %10.5f;  ev = %s'%(x, np.argsort(ev))
            chg_idx = np.argsort(ev)
            tmp0, tmp1 = chg_idx
##            ev0.append(ev[tmp0])
##            ev1.append(ev[tmp1])
            high_state_b0.append( ef[:,tmp0][0]**2 )
            high_state_b1.append( ef[:,tmp0][1]**2 )
            low_state_b0.append( ef[:,tmp1][0]**2 )
            low_state_b1.append( ef[:,tmp1][1]**2 )

            if x <= 50:     # in front-end
                xf_plot.append(x)
                ev0f.append(ev[tmp0])
                ev1f.append(ev[tmp1])
            elif x >= 110:  # in back-end
                xb_plot.append(x)
                ev0b.append(ev[tmp0])
                ev1b.append(ev[tmp1])
            else:           # in middle
                xm_plot.append(x)
                ev0m.append(ev[tmp0])
                ev1m.append(ev[tmp1])
            #endif
        #end of x loop
    elif x_cell == 1:   ## at the begining of 3x3 region
        ####================####
        #### The 3x3 Matrix ####
        H = np.zeros([3,3])
        E0a = 3340
        E0b = 3220
        E0c = 3210
        H[0,1], H[1,0] = [35, 35]
        H[0,2], H[2,0] = [15, 15]

        #############################################
        ## Update the matrix elements, then solve it!
        for x in np.arange( 0, max_x, xinterval ):
            straight_Ea.append(E0a - x)
            straight_Eb.append(E0b + 0.5 * x)
            straight_Ec.append(E0c)
            H[0][0] = straight_Ea[-1]
            H[1][1] = straight_Eb[-1]
            H[2][2] = straight_Ec[-1]
            ev, ef = linalg.eig(H)
            x_plot.append(x)
            chg_idx = np.argsort(ev)
            tmp0 = chg_idx[0]
            tmp1 = chg_idx[1]
            tmp2 = chg_idx[2]
##            ev0.append(ev[tmp0])
##            ev1.append(ev[tmp1])
##            ev2.append(ev[tmp2])
            
            low_state_b0.append( ef[:,tmp0][0]**2 )
            low_state_b1.append( ef[:,tmp0][1]**2 )
            low_state_b2.append( ef[:,tmp0][2]**2 )
            mid_state_b0.append( ef[:,tmp1][0]**2 )
            mid_state_b1.append( ef[:,tmp1][1]**2 )
            mid_state_b2.append( ef[:,tmp1][2]**2 )
            high_state_b0.append( ef[:,tmp2][0]**2 )
            high_state_b1.append( ef[:,tmp2][1]**2 )
            high_state_b2.append( ef[:,tmp2][2]**2 )

            if x <= 50:     # in front-end
                xf_plot.append(x)
                ev0f.append(ev[tmp0])
                ev1f.append(ev[tmp1])
                ev2f.append(ev[tmp2])
            elif x >= 110:  # in back-end
                xb_plot.append(x)
                ev0b.append(ev[tmp0])
                ev1b.append(ev[tmp1])
                ev2b.append(ev[tmp2])
            else:  # in middle
                xm_plot.append(x)
                ev0m.append(ev[tmp0])
                ev1m.append(ev[tmp1])
                ev2m.append(ev[tmp2])
            #endif
        #end of x loop
    #endif
    for y_cell in (0, 1, 2, 3):
        if y_cell == 0:
            ##########################
            ## Plot Eigen-energies
            ax1 = plt.subplot( gs[y_cell, x_cell] )
            ax1 = plt.plot(xf_plot, ev0f, 'b-', xm_plot, ev0m, 'k-', xb_plot, ev0b, 'r-' )
            ax1 = plt.plot(xf_plot, ev1f, 'r-', xm_plot, ev1m, 'k-', xb_plot, ev1b, 'b-' )
            ax1 = plt.plot(x_plot, straight_Ea, 'r:')  # auxiliary line
            ax1 = plt.plot(x_plot, straight_Eb, 'b:')  # auxiliary line
            if x_cell == 1: # 3x3 model
                ax1 = plt.plot(xf_plot, ev2f, 'r-', xm_plot, ev2m, 'k-', xb_plot, ev2b, 'b-' )
                plt.plot(x_plot, straight_Ec, 'k:')  # auxiliary line
            #endif
            xtick_loc_lst = [0, 50, 100, 150]
            ytick_loc_lst = [3150, 3200, 3250, 3300, 3350]
            ax1 = plt.axis(set_energy_axis)
            ax1 = plt.text(30, 3350, 'Eigen-energy', fontsize = title_size) # as title
            ##ax1 = plt.text(40, 3330, 'Eigenstate Energy', fontsize = title_size)
            if right_ylabel:
                #######
                ## Axis-label on the right -- twinx() scheme
                if   x_cell == 0: # 2x2 model
                    plt.xticks(xtick_loc_lst, [])
                    plt.yticks(ytick_loc_lst, [], fontsize = y1tick_size)
                elif x_cell == 1: # 3x3 model
                    plt.yticks(ytick_loc_lst, [], fontsize = y1tick_size)
                    r_yaxis = plt.twinx()
                    r_yaxis.axis(set_energy_axis)
                    r_yaxis.set_yticks(ytick_loc_lst)  # ytick location
                    r_yaxis.set_yticklabels([3150, '', 3250, '', 3350], fontsize = y1tick_size)  # ytick label & font_size
                    plt.xticks(xtick_loc_lst, [])      # this line should be placed under the 'plt.twinx()' line
                    plt.ylabel('Energy (cm-1)', fontsize = ylabel_size )
                #endif
            else:
                #######
                ## Axis-label on the left -- default scheme
                plt.xticks(xtick_loc_lst, [])
                if   x_cell == 0: # 2x2 model
                    plt.yticks(ytick_loc_lst, [3150, '', 3250, '', 3350], fontsize = y1tick_size)
                    plt.ylabel('Energy (cm-1)', fontsize = ylabel_size )
                elif x_cell == 1: # 3x3 model
                    plt.yticks(ytick_loc_lst, [], fontsize = y1tick_size)
                #endif
            #endif
        elif y_cell == 1:
            ##########################
            ## Plot Weights in High-Energy State
            ax2 = plt.subplot( gs[y_cell, x_cell] )
            if   x_cell == 0: # 2x2 model
                plt.plot(x_plot, high_state_b0, 'b-', \
                         x_plot, high_state_b1, 'r-' )
            elif x_cell == 1: # 3x3 model
                plt.plot(x_plot, high_state_b0, 'r-', \
                         x_plot, high_state_b1, 'b-', \
                         x_plot, high_state_b2, 'k-' )
            #endif
            xtick_loc_lst = [0, 50, 100, 150]
            ytick_loc_lst = [0, 0.2, 0.4, 0.6, 0.8, 1]
            plt.axis(set_state_ax)
            plt.text(40, 0.8, 'High-freq.\n   state', fontsize = title_size)
            plt.xticks(xtick_loc_lst, [])
            if right_ylabel:
                #######
                ## Axis-label on the right -- twinx() scheme
                if   x_cell == 0: # 2x2 model
                    ## control the top axis:
                    ## (This part of code should be placed before 'plt.yticks'
                    ##  whenever we need a special control on yticks.)
                    top_ax2 = plt.twiny()
                    top_ax2.axis(set_state_ax)
                    top_ax2.set_xticks([])
                    plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)
                elif x_cell == 1: # 3x3 model
                    ## set y-labels on the right
                    plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)  # original y-axis
                    r_yaxis = plt.twinx()
                    r_yaxis.axis(set_state_ax)
                    r_yaxis.set_yticks(ytick_loc_lst)  # ytick location
                    r_yaxis.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize = y2tick_size)  # ytick label & font_size
                    plt.ylabel('Ratio for each modes', fontsize = ylabel_size )
                    ## control the top axis:
                    ## (Here we don't use 'plt.twiny()', because we already use 'plt.twinx()' )
                    ## move x-axis to the bottom (so there won't be a top x-axis, i.e. no xticks on the top.)
                    ax2.xaxis.set_ticks_position('bottom')
                    plt.xticks(xtick_loc_lst, [])
                #endif
            else:
                #######
                ## Axis-label on the left -- default scheme
                if   x_cell == 0: # 2x2 model
                    plt.yticks(ytick_loc_lst, fontsize = y2tick_size)
                    plt.ylabel('Ratio for each modes', fontsize = ylabel_size )
                    ## control the top axis:
                    ## (This part of code should be placed before 'plt.yticks'
                    ##  whenever we need a special control on yticks.)
                    top_ax2 = plt.twiny()
                    top_ax2.axis(set_state_ax)
                    top_ax2.set_xticks([])
                elif x_cell == 1: # 3x3 model
                    ## control the top axis:
                    ## (This part of code should be placed before 'plt.yticks'
                    ##  whenever we need a special control on yticks.)
                    top_ax2 = plt.twiny()
                    top_ax2.axis(set_state_ax)
                    top_ax2.set_xticks([])
                    plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)
                #endif
            #endif
        elif y_cell == 2 and x_cell == 1:
            ##########################
            ## Plot Weights in Middle-Energy State
            ## (There will be only 3x3 data for "middle-freq state")
            ax3 = plt.subplot( gs[y_cell, x_cell] )
            plt.plot(x_plot, mid_state_b0, 'r-', \
                     x_plot, mid_state_b1, 'b-', \
                     x_plot, mid_state_b2, 'k-' )
            plt.axis(set_state_ax)
            plt.text(20, 0.9, 'Middle-freq. state', fontsize = title_size)

            if right_ylabel:
                #######
                ## Axis-label on the right -- twinx() scheme
                plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)  # original y-axis
                r_yaxis = plt.twinx()
                r_yaxis.axis(set_state_ax)
                r_yaxis.set_yticks(ytick_loc_lst)  # ytick location
                r_yaxis.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize = y2tick_size)  # ytick label & font_size
                plt.ylabel('Ratio for each modes', fontsize = ylabel_size )
                ## control the top axis:
                ## (Here we don't use 'plt.twiny()', because we already use 'plt.twinx()' )
                ## move x-axis to the bottom (so there won't be a top x-axis, i.e. no xticks on the top.)
                plt.xticks(xtick_loc_lst, [])
                ax3.xaxis.set_ticks_position('bottom')
            else:
                #######
                ## Axis-label on the left -- default scheme
                plt.xticks(xtick_loc_lst, [])
                plt.yticks(ytick_loc_lst, fontsize = y2tick_size)
                plt.ylabel('Ratio for each modes', fontsize = ylabel_size)
                ## control the top axis:
                ## (This part of code should be placed before 'plt.yticks'
                ##  whenever we need a special control on yticks.)
                top_ax3 = plt.twiny()
                top_ax3.axis(set_state_ax)
                top_ax3.set_xticks([])
            #endif
        elif y_cell == 3:
            ##########################
            ## Plot Weights in Low-Energy State
            ax4 = plt.subplot( gs[y_cell, x_cell] )
            if   x_cell == 0: # 2x2 model
                plt.plot(x_plot, low_state_b0, 'b-', \
                         x_plot, low_state_b1, 'r-' )
                plt.text(43, 0.8, 'Low-freq.\n   state', fontsize = title_size)
            elif x_cell == 1: # 3x3 model
                plt.plot(x_plot, low_state_b0, 'r-', \
                         x_plot, low_state_b1, 'b-', \
                         x_plot, low_state_b2, 'k-' )
                plt.text(28, 0.9, 'Low-freq. state', fontsize = title_size)
            #endif
            plt.axis(set_state_ax)

            if right_ylabel:
                #######
                ## Axis-label on the right -- twinx() scheme
                plt.xticks(xtick_loc_lst, fontsize = xtick_size)
                plt.xlabel('x', fontsize = ylabel_size)
                ## control the top axis:
                ## (Here we don't use 'plt.twiny()', because we already use 'plt.twinx()' )
                ## move x-axis to the bottom (so there won't be a top x-axis, i.e. no xticks on the top.)
                ax4.xaxis.set_ticks_position('bottom')
                if   x_cell == 0: # 2x2 model
                    plt.xticks(xtick_loc_lst, fontsize = xtick_size)
                    plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)
                elif x_cell == 1: # 3x3 model
                    plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)  # original y-axis
                    r_yaxis = plt.twinx()
                    r_yaxis.axis(set_state_ax)
                    r_yaxis.set_yticks(ytick_loc_lst)  # ytick location
                    r_yaxis.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1], fontsize = y2tick_size)  # ytick label & font_size
                    plt.ylabel('Ratio for each modes', fontsize = ylabel_size )
                #endif
            else:
                #######
                ## Axis-label on the left -- default scheme
                plt.xticks(xtick_loc_lst, fontsize = xtick_size)
                plt.xlabel('x', fontsize = ylabel_size)
                if   x_cell == 0: # 2x2 model
                    plt.yticks(fontsize = y2tick_size)
                    plt.ylabel('Ratio for each modes', fontsize = ylabel_size)
                    ## control the top axis:
                    ## (This part of code should be placed before 'plt.yticks'
                    ##  whenever we need a special control on yticks.)
                    top_ax4 = plt.twiny()
                    top_ax4.axis(set_state_ax)
                    top_ax4.set_xticks([])
                elif x_cell == 1: # 3x3 model
                    ## control the top axis:
                    ## (This part of code should be placed before 'plt.yticks'
                    ##  whenever we need a special control on yticks.)
                    top_ax4 = plt.twiny()
                    top_ax4.axis(set_state_ax)
                    top_ax4.set_xticks([])
                    plt.yticks(ytick_loc_lst, [], fontsize = y2tick_size)
                #endif
            #endif



###############################################
###############################################
##
## Tip 1. To eliminate the xticks of upper axis.
##   We have two ways:
##   (1) move the x-axis
''' [Code Example]
         ax4 = plt.subplot( gs[y_cell, x_cell] )
         plt.plot([1,2,3], [4,5,6])
         plt.xticks(xtick_loc_lst, fontsize = xtick_size)
         plt.xlabel('x', fontsize = ylabel_size)
         ax4.xaxis.set_ticks_position('bottom')
'''
##   (2) using 'plt.twiny()' to access the top xticks
''' [Code Example]
         ax3 = plt.subplot( gs[y_cell, x_cell] )
         plt.plot([1,2,3], [4,5,6])
         ## control the top axis:
         ## (The following part of code should be placed before 'plt.yticks'
         ##  whenever we need a special control on yticks.)
         top_ax3 = plt.twiny()
         top_ax3.axis(set_state_ax)
         top_ax3.set_xticks([])
'''
##
## Tip 2. To 
###############################################
###############################################



####=============================####
#### Control the Border of Plots ####
'''
We can consider the Figure as a peice of paper.
The 'gridspec' is the rectangular grid on that paper.
To control the location of 'gridspec', we need to introduce two kind of quantities:
  border & spacing.

The border decides which locations (on the figure) should 'gridspec' bigins and ends.
For the horizontal border, we have 'left' and 'right' parameters.
For the vertical border, we have 'top' and 'bottom' parameters.

There is another type of border which is between each subplots.
And we called it spacing.
For the horizontal and vertical spacing, we have 'wspace' and 'hspace' respectively.
'''

#########################
## Set horizontal border: 
if right_ylabel:
    gs.update( left=0.02, right=0.85 )
else:
    gs.update( left=0.15, right=0.98 )
#endif
## Placing a horizontal axis for this figure.
## Total width is 1, the figure's left-edge is 0.0, right-edge is 1.0
## Parameter 'left'  indicates the horizontal location where 'gridspec' begins,
##           'right' indicates the horizontal location where 'gridspec' ends.

#########################
## Set vertical border: 
gs.update( bottom=0.05, top=0.95 )
## Placing a vertical axis for this figure.
## Total width is 1, the figure's bottom-edge is 0.0, top-edge is 1.0
## Parameter 'bottom'  indicates the vertical location where 'gridspec' begins,
##           'top' indicates the vertica location where 'gridspec' ends.

#########################
## Set spacing: 
gs.update( wspace = 0.15 , hspace = 0.15 )
## Parameter 'wspace' is the horizontal space between each subplot,
##           'hspace' is the vertical space between each subplot.


####============================####
#### Control the Size of Figure ####
cm2inch = 0.393700787
ratio_w_to_h = 4.0 / 2.3
fig_width  = 8.5 * cm2inch  #fig.get_size_inches()[0]
fig_height = fig_width * ratio_w_to_h #fig.get_size_inches()[1]

fig.set_figheight(fig_height)  #unit: inch
fig.set_figwidth(fig_width)   #unit: inch

plt.savefig('Plot_2x2_and_3x3_toy_model.png', dpi=600, format='png')
##plt.show()
print "Figure complete."

