import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
# Set plotting parameters

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams['ps.fonttype'] = 42

pixels_needed = 500 / 0.65
fontprops = fm.FontProperties(size=18)

def gaussian_heat_plot(x,
                       y,
                       vmin_val,
                       vmax_val,
                       title_text,
                       plot_contour=False,
                       title=False,
                       save_name='',
                       save=False):
    """Gaussian smoothed plots for gene expression on puck"""

    data = np.array([np.array([i, j]) for i, j in zip(x, y)])
    k = kde(data.T)

    xi, yi = np.mgrid[0:6000:nbins * 1j, 0:6000:nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    print("real life size:", (min(x) - max(x) / nbins * 0.65))

    Z, xedges, yedges = np.histogram2d(x, y, bins=nbins)

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim([0, 6000])
    ax.set_ylim([0, 6000])
    if title:
        ax.set_title(title_text)
    Z = gaussian_filter(Z, sigma=0.5)
    this_plot = plt.pcolormesh(
        xedges[:-1],
        yedges[:-1],
        Z.T,
        cmap=reversed_color_map,
        shading="gouraud",
        vmin=vmin_val,
        vmax=vmax_val,
    )
    if plot_contour:
        ax.contour(xi,
                   yi,
                   zi.reshape(xi.shape),
                   levels=15,
                   linewidths=0.5,
                   colors="white")
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "right",
        size="7%",
        pad=0.2,
    )
    cbar = fig.colorbar(this_plot, cax=cax, ticks=[])
    cbar.ax.tick_params(labelsize=20)
    ax.set_axis_off()

    cbar.ax.text(0.5,
                 -0.01,
                 str(vmin_val),
                 transform=cbar.ax.transAxes,
                 va="top",
                 ha="center")
    cbar.ax.text(0.5,
                 1.0,
                 str(vmax_val),
                 transform=cbar.ax.transAxes,
                 va="bottom",
                 ha="center")
    cbar.outline.set_visible(False)
    scalebar = AnchoredSizeBar(
        ax.transData,
        pixels_needed,
        "",
        "lower right",
        pad=0.1,
        color="black",
        frameon=False,
        size_vertical=1,
        fontproperties=fontprops,
    )

    ax.add_artist(scalebar)

    fig.canvas.draw()

    if save == True:
        fig.savefig(f"./plots/{title_text}_{save_name}.pdf",
                    bbox_inches='tight')

        
# Plot cortex/medulla regions
def plot_cortex_medulla(adata, save_name = False):
    fig,ax = plt.subplots(figsize=(5,5))
    x = adata.obs[adata.obs['cortex_medulla'] == 'cortex'].x
    y = adata.obs[adata.obs['cortex_medulla'] == 'cortex'].y
    plt.scatter(x,y,label='cortex',color="#0F4C81",s=10)

    x = adata.obs[adata.obs['cortex_medulla'] == 'medulla'].x
    y = adata.obs[adata.obs['cortex_medulla'] == 'medulla'].y
    plt.scatter(x,y,label='medulla',color="#FFD470",s=10)
    
    scalebar = AnchoredSizeBar(
        ax.transData,
        pixels_needed,
        "",
        "lower right",
        pad=0.1,
        color="black",
        frameon=False,
        size_vertical=1,
        fontproperties=fontprops,
    )

    ax.add_artist(scalebar)
    plt.xlim([0,5000])
    plt.ylim([0,5000])
    
    if save_name != False:
        plt.savefig(f'{save_name}.pdf')
    
        