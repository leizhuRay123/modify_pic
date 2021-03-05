from const import *
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

dirc = sys.argv[1]

def plot_field(dirc,i):
    info = sr.read_const(dirc)
    ffs = sr.Get_file(prefix = '',dirc = dirc)
    extent = list(np.array([info['xmin'],info['xmax'],info['ymin'],info['ymax']])/info['de0'])
    normed_ne = 1e-1*info['n0']
    a = sdf.read(ffs[i])
    xylims = [extent[0]+40,-40,extent[2]/2,extent[3]/2]
    fig,axs = df.Create_Figure2(w=2,h=3,sp_h=-1)
#     fig.suptitle(str(a.Header['time']/info['T0']))
    mvp,xx,yy = sr.get_mvp(sdffile=ffs[i],avg=False)
    
    ax = axs[0]
    df.draw_field_snapshot(ax = ax, data = a.Magnetic_Field_Bx.data.T,extent=extent)
    ax.contour(xx/info['de0'],yy/info['de0'],mvp.T,levels = 40)
    df.Axis_set(ax, axesname=['','','Bx ' + str(a.Header['time']/info['T0'])[0:4]+'$T_0$'],yticklabal=[],xylims = xylims)
    ax.set_xticks([])
    ax.set_yticks([])
    ax = axs[1]
    # jz = a.Current_Jz
    df.draw_field_snapshot(ax = ax, data = a.Magnetic_Field_By.data.T,cmap='bwr',extent=extent)
    df.Axis_set(ax,axesname=['','','By'],xylims = xylims)
    # axis_xticks = ax.get_xticks()
    ax.set_xticks([])
    ax.set_yticks([])


    ax = axs[2]
    df.draw_field_snapshot(ax = ax, data = a.Magnetic_Field_Bz.data.T/info['j0']/info['de0'],extent=extent)
    df.Axis_set(ax,axesname=['','','Bz'],xylims = xylims)
    ax.set_xticks([])
    ax.set_yticks([])
    # ax = axs[3]
    # df.draw_field_snapshot(ax = ax, data = a.Derived_Number_Density.data.T/(info['n0']),extent=extent)
    # df.Axis_set(ax,axesname=['','','ne'])

    ax = axs[3]
    # jx = a.Current_Jx
    df.draw_field_snapshot(ax = ax, data = a.Electric_Field_Ez.data.T/c,cmap='bwr',extent=extent)
    df.Axis_set(ax,axesname=['','','Ez'],xylims = xylims)
    ax.set_xticks([])
    ax.set_yticks([])

    ax = axs[4]
    # jx = a.Current_Jx
    ne = sr.Get_field_variable(a,var='Number_Density').data

    df.draw_field_snapshot(ax = ax, data = np.log10(ne.T/(normed_ne)),\
                           cmap='jet',extent=extent,sym=False,caxis = [-1,1.5])

    df.Axis_set(ax,axesname=['','','ne'],xylims = xylims)
    ax.set_xticks([])
    ax.set_yticks([])

    ax = axs[5]
    # jx = a.Current_Jx
    varname = 'Current_Jz0'
    jz = sr.Get_field_variable(a,varname).data
#     nx,ny = jz.shape
#     print(jz)
#     ax.plot(jz[:,ny//2])
    df.draw_field_snapshot(ax = ax, data = jz.T,\
                           cmap='bwr',extent=extent,caxis = [-1,1])

    df.Axis_set(ax,axesname=['','',varname],xylims = xylims)
    ax.set_xticks([])
    ax.set_yticks([])
    
    
    
for i in range(len(ffs)):
    if (rank == i % size):
        plot_field(dirc,i)
        plt.savefig(dirc + 'FE'+str(i)+'.png',dpi = 200)
        plt.close('all')
    else:
        continue;
