# Output reading test

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap



# Colormap
# =====================
cdict1 = {'red':  ((0.0 , 0.0, 0.0),
                   (0.5, 0.0, 0.0),
                   (0.75 , 1.0, 1.0),
                   (1.0 , 1.0, 1.0)),

         'green': ((0.0 , 1.0, 1.0),
                   (0.25, 0.0, 0.0),
                   (0.75, 0.0, 0.0),
                   (1.0 , 1.0, 1.0)),

         'blue':  ((0.0 , 1.0, 1.0),
                   (0.25 , 1.0, 1.0),
                   (0.5 , 0.0, 0.0),
                   (1.0 , 0.0, 0.0))
        }


CMAP = LinearSegmentedColormap('StokesFD', cdict1)
plt.register_cmap(cmap=CMAP)
plt.set_cmap('StokesFD')
#plt.set_cmap("gray")


plt.pcolor([[0,1],[0,1]],[[0,1],[0,1]],[[0,1],[0,1]],vmin=-1, vmax=2)


plt.colorbar(orientation="horizontal")
plt.show()
        
        
        
        
        
        
        #plt.savefig("HFac" + str(round(thisHFac)) + "_C" + str(round(thisCohesion)) + ".png", dpi=500)
        
        













