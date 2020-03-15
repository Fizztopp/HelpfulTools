import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl


def test_lin(m,x):
    return m*x
    
def test_quad(m,x):
    return m*x**2

aa=1.5    
    
II = np.linspace(1,100,100)

fig2 = plt.figure(2)
gs2 = gridspec.GridSpec(1, 1)

X=np.log(II)
Y_LIN=np.log(test_lin(aa,II)) 
Y_QUAD=np.log(test_quad(aa,II)) 

fit_lin = np.polyfit(X, Y_LIN, 1)
p_lin = np.poly1d(fit_lin)
fit_lin = lambda x: np.exp(p_lin(np.log(x)))
pdx_lin = np.polyder(p_lin,1)
print("linearized function test_lin(): fit(test_lin(m,x))="+str(p_lin)+" --> linear")

fit_quad = np.polyfit(X, Y_QUAD, 1)
p_quad = np.poly1d(fit_quad)
fit_quad = lambda x: np.exp(p_quad(np.log(x)))
pdx_quad = np.polyder(p_lin,1)
print("linearized function test_quad(): fit(test_quad(m,x))="+str(p_quad)+" --> quadratic")

ax11 = fig2.add_subplot(gs2[0,0])
ax11.set_xlabel(r'$\mathrm{x}$',fontsize=24)
ax11.set_ylabel(r'$\mathrm{f(x)}$',fontsize=24)

ax11.loglog(II,test_lin(aa,II), "bx",label="test_lin(m,x)")
ax11.loglog(II,test_quad(aa,II), "rx",label="test_quad(m,x)")
#ax11.axvline(x=10.0, ymin=1.0, ymax = 10, c="k", linestyle="--")
#ax11.axvline(x=100.0, ymin=1.0, ymax = 10, c="k", linestyle="--")
ax11.loglog(II,fit_lin(II), "b",label="fit(test_lin(m,x))")
ax11.loglog(II,fit_quad(II), "r",label="fit(test_quad(m,x))")

#ax11.set_xlim(0.9,12)
#ax11.set_ylim(1.0,100)  
plt.legend()
 
plt.tight_layout()
plt.show()    