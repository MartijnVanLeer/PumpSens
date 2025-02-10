#%%
import matplotlib.pyplot as plt

def profile(ax, confined):
# sky
# 
    if confined:
        sky = plt.Rectangle((-25, 0), width=250, height=5, fc="black", zorder=0, alpha=1)
    else:
        sky = plt.Rectangle((-25, 0), width=250, height=5, fc=np.array([100, 100, 100]) / 255, zorder=0, alpha=0.7)

    bot = plt.Rectangle((-25, -35), width=250, height=5, fc="black", zorder=0, alpha=1)


    # Aquifer:
    ground = plt.Rectangle(
        (-25, -10),
        width=250,
        height=10,
        fc=np.array([209, 179, 127]) / 255,
        zorder=0,
        alpha=0.9,
    )


    # Aquifer 2:
    ground2 = plt.Rectangle(
        (-25, -30),
        width=250,
        height=10,
        fc=np.array([209, 179, 127]) / 255,
        zorder=0,
        alpha=0.9,
    )

    # Confining bed:
    confining_unit = plt.Rectangle(
        (-25, -20),
        width=250,
        height=10,
        fc=np.array([100, 100, 100]) / 255,
        zorder=0,
        alpha=0.7,
    )
    well = plt.Rectangle(
        (-1.5, -30), width=3, height=35, fc=np.array([200, 200, 200]) / 255, zorder=1
    )

    # Screen for the well:
    screen = plt.Rectangle(
        (-1.5, -30),
        width=3,
        height=10,
        fc=np.array([200, 200, 200]) / 255,
        alpha=1,
        zorder=2,
        ec="k",
        ls="--",
    )
    screen.set_linewidth(2)


    obswell = plt.Rectangle(
        (198.5, -30), width=3, height=35, fc=np.array([200, 200, 200]) / 255, zorder=1
    )

    obsscreen = plt.Rectangle(
        (198.5, -30),
        width=3,
        height=10,
        fc=np.array([200, 200, 200]) / 255,
        alpha=1,
        zorder=2,
        ec="k",
        ls="--",
    )
    obsscreen.set_linewidth(2)


    obsscreen2 = plt.Rectangle(
        (198.5, -10),
        width=3,
        height=10,
        fc=np.array([200, 200, 200]) / 255,
        alpha=1,
        zorder=2,
        ec="k",
        ls="--",
    )
    obsscreen2.set_linewidth(2)

    ax.add_patch(obsscreen2)
    ax.add_patch(obsscreen)
    ax.add_patch(obswell)
    ax.add_patch(screen)
    ax.add_patch(well)
    ax.add_patch(confining_unit)
    ax.add_patch(ground2)
    ax.add_patch(ground)
    ax.add_patch(bot)
    ax.add_patch(sky)

    ax.set_yticks([0,-10,-20,-30 ])
    ax.set_xlim([-25, 225])
    ax.set_ylim([-35, 5])
    return ax

def add_text(ax,confined):
    x = 87
    if confined:
        ax.text(x=x, y=2, s=r"No flow", c = 'white')
    else: 
        ax.text(x=x-20, y=2, s=r"Semiconfining layer", c = 'black')

    ax.text(x=x, y=-33, s=r"No flow", c = 'white')
    ax.text(x=x, y=-26, s=r"Aquifer 2")
    ax.text(x=x, y=-16, s=r"Aquitard")
    ax.text(x=x, y=-6, s=r"Aquifer 1")
    ax.text(x=-20, y=7, s=r"Extraction well")
    ax.text(x=160, y=7, s=r"Observation wells")

    ax.set_xlabel("Distance [m]")
    ax.set_ylabel("Relative height [m]")
    return ax

def get_color(width):
    if width == 10:
        color = '#8a0000'
    elif width == 6:
        color = 'red'
    elif width == 3:
        color = '#ff8080'
    return color

def arrows_horizontal(ax, layer, width, linestyle = 'solid'):
    z = -25 if layer == 2 else -5
    fc = 'none' if linestyle == 'dashed' else get_color(width)
    arrow2 = plt.Arrow(30,z,140,0, width = width, color = get_color(width), ec =  'black', linestyle = linestyle, fc = fc)
    ax.add_patch(arrow2)
    return ax


def arrows_vertical(ax, well, width):
    fac = 6.25
    x = 170 if well == 'obs' else 30

    arrow1 = plt.Arrow(x,-23,0,16, width = width*fac, color = get_color(width), ec =  'black')
    ax.add_patch(arrow1)
    return ax

def arrows_sf(ax,well, width):
    fac = 6.25
    x = 170 if well == 'obs' else 30
    ly = -4 if well == 'obs' else -6
    arrow1 = plt.Arrow(x,3,0,ly, width = width*fac, color = get_color(width), ec =  'black')
    arrow2 = plt.Arrow(100,3,0,-5, width = width*fac, color = get_color(width), ec =  'black')
    ax.add_patch(arrow1)
    ax.add_patch(arrow2)
    return ax



def plot():
    fig, axs = plt.subplots(3,4, constrained_layout = True, sharex = True, sharey = True)
    fig.set_size_inches(17.15/2.54,5)
    for col in [0,1]:
        for row in range(3):
            axs[row,col] = profile(axs[row,col], True)
            axs[row,col+2] = profile(axs[row,col+2], False)
            axs[row,col].set_yticks([])
            axs[row,col].set_xticks([])

    wmax = 10
    wmid = 6
    wmin = 3
    #confined
    axs[0,0].set_ylabel('K1 < K2')
    axs[1,0].set_ylabel('K1 = K2')
    axs[2,0].set_ylabel('K1 > K2')
    axs[0,2].set_ylabel(' ')
    axs[1,2].set_ylabel(' ')
    axs[2,2].set_ylabel(' ')
    for col in [0,2]:
        axs[0,col].set_title('Early')
        axs[0,col+1].set_title('Late')
        #early K1 < K2
        arrows_horizontal(axs[0,col], layer =2, width = wmax)
        arrows_vertical(axs[0,col], well = 'obs', width = wmax)
        #early K1 = K2
        arrows_horizontal(axs[1,col], layer =2, width = wmid)
        arrows_vertical(axs[1,col], well = 'obs', width = wmid)
        arrows_horizontal(axs[1,col], layer =0, width = wmid)
        arrows_vertical(axs[1,col], well = 'pump', width = wmid)
        #early K1 > K2
        arrows_horizontal(axs[2,col], layer =0, width = wmax)
        arrows_vertical(axs[2,col], well = 'pump', width = wmax)

    #confined
    #lateK1 < K2
    arrows_horizontal(axs[0,1], layer =2, width = wmax)
    arrows_vertical(axs[0,1], well = 'obs', width = wmax)
    arrows_horizontal(axs[0,1], layer =0, width = wmin)
    arrows_vertical(axs[0,1], well = 'pump', width = wmin)
    #late K1 = K2
    arrows_horizontal(axs[1,1], layer =2, width = wmid)
    arrows_vertical(axs[1,1], well = 'obs', width = wmid)
    arrows_horizontal(axs[1,1], layer =0, width = wmid)
    arrows_vertical(axs[1,1], well = 'pump', width = wmid)
    #late K1 > K2
    arrows_horizontal(axs[2,1], layer =0, width = wmax)
    arrows_vertical(axs[2,1], well = 'pump', width = wmax)
    arrows_horizontal(axs[2,1], layer =2, width = wmin)
    arrows_vertical(axs[2,1], well = 'obs', width = wmin)

    #semiconfined
    #late K1 <K2 zelfde als early
    arrows_horizontal(axs[0,3], layer =2, width = wmax)
    arrows_vertical(axs[0,3], well = 'obs', width = wmax)
    arrows_horizontal(axs[0,3], layer =0, width = wmin, linestyle = 'dashed')
    arrows_vertical(axs[0,3], well = 'pump', width = wmin)
    #K1 = K2
    arrows_horizontal(axs[1,3], layer =2, width = wmid)
    arrows_vertical(axs[1,3], well = 'obs', width = wmid)
    arrows_horizontal(axs[1,3], layer =0, width = wmin)
    arrows_vertical(axs[1,3], well = 'pump', width = wmid)

    arrows_horizontal(axs[2,3], layer =2, width = wmid)
    arrows_vertical(axs[2,3], well = 'obs', width = wmid)
    arrows_horizontal(axs[2,3], layer =0, width = wmid)
    arrows_vertical(axs[2,3], well = 'pump', width = wmax)

    for row in range(3):
        arrows_sf(axs[row,3],'pump', width = wmin)
        arrows_sf(axs[row,3],'obs', width = wmin)
            
    
    return fig, axs

fig, axs = plot()