# %%
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams['pgf.preamble'] = [r'\usepackage{url}', ]
from matplotlib import pyplot as plt
import numpy as np
from pandas import read_csv
import webbrowser
import pandas as pd

plot_points_without_errors = True
color_quantity = "Year"
colormap = "viridis"
marker_by_type = True
mmin = 0.01

if marker_by_type:
    markerdict = {"OB Association": "d", "Young Cluster": "o", "MW Field": ">", "MW Bulge": "v","Globular Cluster": "X", "MW Nuclear Cluster": "^"}

data = read_csv("alphaplot.csv",skip_blank_lines=True)
data = data[~np.isnan(data["Slope (Salpeter=2.35)"])]

references = data["Reference"]
system = data["System"]
Z = data["Metallicity [Z/H]"]
#cut = (system=="ONC") #(Z<=-0.5) #np.ones(len(data),dtype=np.bool) #(system=="ONC")#*(year <= 2003)
#cut = np.ones(len(data))
#data = data[cut]
references = data["Reference"]
types = data["Class"]
markers = np.array([markerdict[t] for t in types])
year = []
for r in references:
    try: 
        year.append(int(r[:4]))
    except:
        year.append(1995)
year = np.array(year)
#year = np.array([int(r[:4]) if type(r)==str else -1 for r in references])
slope = np.array(-data["Slope (Salpeter=2.35)"]+1)
slope += np.random.normal(size=(len(slope),))*0.02
mlow = np.array(data["Lower mass (Msun)"])
mhi = np.array(data["Upper mass (Msun)"])
mmed = np.array((mlow*mhi)**0.5)
slope_error = np.array(data["Slope uncertainty"])
slope_upper = np.array(-data["Upper limit (1 sigma)"]+1)
slope_lower = np.array(-data["Lower limit (1 sigma)"]+1)

# %%

fig, ax = plt.subplots(1,1,figsize=(6,3))

if color_quantity=="Year":
    colors = plt.get_cmap(colormap)((year-year[year>0].min()).clip(0,1e100)/(year-year[year>0].min()).max())
elif color_quantity=="Metallicity":
    colors = plt.get_cmap(colormap)((Z-Z.min())/(Z.max()-Z.min()))
else:
    colors=np.zeros((len(data),3))


ebar_lw = 0.3 # linewidth of errorbars
# first plot asymmetric errors where available
i = np.isfinite(slope_lower) * np.isfinite(slope_upper)
ax.errorbar(mmed[i],slope[i],xerr=np.array([mmed[i]-mlow[i],mhi[i]-mmed[i]]), yerr=np.array([np.abs(slope-slope_lower)[i],np.abs(slope_upper-slope)[i]]), ls='',capsize=0,lw=ebar_lw,c='black',marker=None,ecolor='grey')
for m in np.unique(markers):
    ax.scatter(mmed[i*(markers==m)],slope[i*(markers==m)],c=colors[i*(markers==m)],s=10,zorder=20,marker=m,lw=0) #xerr=[mmed[i]-mlow[i],mhi[i]-mmed[i]], yerr=[(slope-slope_lower)[i],(slope_upper-slope)[i]], ls='',capsize=0,lw=0.5,c='black',marker='o',markersize=1)

# now symmetric errors
i = np.isfinite(slope_error)# * np.isfinite(slope_upper)
ax.errorbar(mmed[i],slope[i],xerr=np.array([mmed[i]-mlow[i],mhi[i]-mmed[i]]), yerr=slope_error[i], ls='',capsize=0,lw=ebar_lw,c='black',marker=None,markersize=0,ecolor='grey')
for m in np.unique(markers):
    ax.scatter(mmed[i*(markers==m)],slope[i*(markers==m)],c=colors[i*(markers==m)],s=10,zorder=10,marker=m,lw=0)

# and now ones without any errorbars :(
if plot_points_without_errors:
    i = np.isnan(slope_error)
    ax.errorbar(mmed[i],slope[i],xerr=np.array([mmed[i]-mlow[i],mhi[i]-mmed[i]]), ls='',capsize=0,lw=ebar_lw,color='black',marker='o',markersize=0,ecolor='grey')
    
    for m in np.unique(markers):
        ax.scatter(mmed[i*(markers==m)],slope[i*(markers==m)],c=colors[i*(markers==m)],s=10,zorder=2,marker=m,lw=0)


#colorbar
if color_quantity=="Year":
    import matplotlib.ticker as ticker
    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'%d'%int(x)

    sc = ax.scatter(np.zeros(len(mmed)+1),-100*np.ones(len(mmed)+1),c=np.int_(list(year)+[year[year>0].min()]),s=10,cmap=colormap)
    plt.colorbar(sc,format=ticker.FuncFormatter(fmt),pad=0,label="Year")
elif color_quantity=="Metallicity":
    sc = ax.scatter(np.zeros(len(mmed)),-100*np.ones(len(mmed)),c=Z,s=10,cmap=colormap)
    plt.colorbar(sc,label="Z",pad=0)


if len(np.unique(data["System"].values))==1: 
    ax.set_title(data["System"][0])

# plot dummy points for URLS
# df = pd.DataFrame({'x': mmed,
#                    'y': slope,
#                    'link': [r"https://ui.adsabs.harvard.edu/abs/"+r for r in references]})

# def on_pick(event):
#     url = df.link.iloc[event.ind[0]]
#     webbrowser.open_new_tab(url)

# ax.scatter(x=df.x, y=df.y,picker=5,s=1,alpha=0)

urls = [r"https://ui.adsabs.harvard.edu/abs/"+r.split(";")[0] for r in references]
# import matplotlib.patches as patches
sc = ax.scatter(mmed,slope,alpha=1e-6,zorder=10000)
sc.set_urls(urls)
for x,y,u in zip(mmed,slope,urls):
    txt = ax.text(x, y, "o", url=u, alpha=1e-6, fontsize=3,bbox=dict(boxstyle='circle',url=u,alpha=1e-6),ha='center',va='center',zorder=10000) #bbox = dict(color='w', alpha=0.01, url=u)
    #ax.add_patch(patches.Rectangle((x,y),0.1,0.1,url=u))#, url=u, alpha=1e-6, fontsize=3,bbox=dict(boxstyle='circle',url=u,alpha=1e-6),ha='center',va='center',zorder=10000)

peak_pos = 0.2
ax.plot([mmin,300],[0,0],color='black',zorder=-1000,lw=0.5,ls='solid')
ax.fill_between([0.1,0.3],[-0.3,-0.3],[0.3,0.3],color='purple',zorder=-1000,alpha=0.3,label="Typical Peak/Plateau",linewidth=0)
#ax.plot([0.2,0.2],[-10,10],ls='dashed',color='grey',zorder=-1000,lw=0.5)
#ax.text(0.0045,-0.2,'Peak/Plateau',color='black',fontsize=6)
#ax.text(0.0045,0.4,'$\sigma$',color='black',fontsize=6)
#ax.text(0.0045,-0.6,'$\sigma$',color='black',fontsize=6)
#ax.arrow(0.004,0,0,1,width=0.0001,head_length=0.1,head_width=0.001,length_includes_head=True,facecolor='black',lw=0)
#ax.scatter(0.004,0,color='black',s=4)
#ax.arrow(0.004,0,0,-1,width=0.0001,head_length=0.1,head_width=0.001,length_includes_head=True,facecolor='black',lw=0)

ax.plot([mmin,300],[-1.35,-1.35],ls='-.',color='black',zorder=-1000,lw=0.5)
ax.text(mmin*2.2,-1.55,"Salpeter Slope (-1.35)",color='black',fontsize=6)
# ax.plot([mmin,300],[-1.,-1.],ls='dotted',color='grey',zorder=-1000,lw=0.5,label="Log-normal $\pm \sigma$")
# ax.plot([mmin,300],[1.,1.],ls='dotted',color='grey',zorder=-1000,lw=0.5)
# ax.plot([0.099,0.099],[-10.,10.],ls='dotted',color='grey',zorder=-1000,lw=0.5)
# ax.plot([0.7096,0.7096],[-10.,10.],ls='dotted',color='grey',zorder=-1000,lw=0.5)

mgrid = np.logspace(-3,3,1000)

def analytic_imf_slope(mgrid,model="Chabrier (2005)"):
    slope = np.zeros_like(mgrid)
    if "Chabrier" in model:
        if "2005" in model: 
            logmc = np.log10(0.2) # 
            sigma = 0.55
        else:
            logmc = np.log10(0.08)
            sigma = 0.69
        slope = -(np.log10(mgrid)-logmc)/(sigma**2) / np.e
        slope[mgrid>1] = -1.35
    elif "Kroupa" in model:
        slope[mgrid<0.08] = 0.7
        if "2002" in model:
            slope[(mgrid>0.08)*(mgrid<0.5)] = -0.3
            slope[mgrid>0.5] = -1.3
        elif "2001" in model:
            slope[(mgrid>0.08)*(mgrid<0.5)] = -0.8
            slope[(mgrid>0.5)*(mgrid<1)] = -1.7
            slope[mgrid>1] = -1.3
    
    slope[mgrid>120] = np.nan
    slope[mgrid<0.01] = np.nan

    return slope

modelcolors = "red", "darkblue"
models = "Kroupa (2002)","Chabrier (2005)"
for c, model in zip(modelcolors,models):
    slope_model = analytic_imf_slope(mgrid,model)
    ax.plot(mgrid,slope_model,zorder=-1000,label=model,color=c)


for t, m in markerdict.items():
    plt.scatter([-10,-10],[-10,-10],marker=m,lw=0,color='black',label=t,s=10)

ax.set(xscale='log',xlabel=r"Stellar Mass ($M_\odot$)", ylabel=r"IMF Slope $\Gamma_{\rm IMF}$",xlim=[mmin,300],ylim=[-3,3.5])
ax.legend(fontsize=6,labelspacing=0)
#fig.canvas.mpl_connect('pick_event', on_pick)
plt.savefig("IMF_AlphaPlot.pdf",bbox_inches='tight')


