# -*- coding: utf-8 -*-

#  GalerkinAppr 0.7
#  
#  Copyright 2020: Juan Carlos del Caño
#                  Prof. at the "Escuela de Ingenierias Industriales"
#                  University of Valladolid (Spain)
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# version 0.7 - Primera version con GUI.


# Importacion de librerias 
import numpy as np
from numpy import exp, sin, cos, log  # para ponerlo asi en las ff
import sympy as sp
from tkinter import *
from tkinter import ttk, messagebox, filedialog

import matplotlib
from matplotlib import pyplot as plt
matplotlib.use('TkAgg')  # mac requiere hacer esto explicitamente
from matplotlib import gridspec
from os import path, kill, getppid, _exit
import signal



class Nodo:
    def __init__(self, x=0., y=0., ux=0., uy=0., elemn=[], sigma=[0.,0.,0.]):
        # solo almacena info, no se necesitaria mas (?)
        self.x, self.y = x,y
        self.ux, self.uy = ux, uy
        self.elemn = elemn
        self.sigma = sigma

class Carga:
    def __init__(self, contorno=[], puntual=[], volumen=[], termica=[]):
        # En principio habra una sola instancia. Cada lista tiene en [0]
        # el numero de cargas de ese tipo y luego sub listas con las cargas.
        # Para contorno [nt ó xy, nodo, nodo, nodo, p1, p1, p1, p2, p2, p2]
        # Para puntual [nodo, Fx, Fy]
        # Para volumen [string de la funcion Xx, idem Xy]
        # Para termica [string de la funcion T]
        
        self.contorno= contorno
        self.puntual = puntual
        self.volumen = volumen
        self.termica = termica
    

class Elem:
    # Guardo como atributos de clase los coeficientes de las ff geometricas (H)
    # las coord normalizadas de los nodos (-1,1 o similares) y
    # las coord normalizadas de los puntos de Gauss y sus pesos.
    # No quiero que se instancien ya que son igual para cada tipo de elemento.
    
    
    # Coef de las ff geometricas H. Hay que multiplicar asi:
    # [1, x_, y_, x_**2, x_*y_, etc] *[coef[][]] = [H0(x_y_), H1(x_y_), etc ] 
    # para obtener los valores de las ff H en (x_,y_) ¡que son coor normalizadas!
    coefH=[None, None]
    
    # en [0] lo del cuadrado_8
    a,b= 0.25, 0.5
    coefH[0]=np.array([
    [-a, b,-a, b,-a, b,-a, b],
    [ 0, 0, 0, b, 0, 0, 0,-b],
    [ 0,-b, 0, 0, 0, b, 0, 0],
    [ a,-b, a, 0, a,-b, a, 0],
    [ a, 0,-a, 0, a, 0,-a, 0],
    [ a, 0, a,-b, a, 0, a,-b],
    [-a, b,-a, 0, a,-b, a, 0],
    [-a, 0, a,-b, a, 0,-a, b]   ])

    # en [1] lo del triangulo_6
    r3=np.sqrt(3.)
    coefH[1]=np.array([
    [ 1.00,    0.00,    0.00,    0.00,    0.00,   0.00 ],
    [ 0.00,    2.00,   -0.50,    0.00,    0.50,  -2.00 ],
    [  -r3,    2/r3,  -0.5/r3,   0.00, -0.5/r3,   2/r3 ],
    [ 0.00,    0.00,    0.50,   -1.00,    0.50,   0.00 ],
    [ 0.00,   -2/r3,    1/r3,    0.00,   -1/r3,   2/r3 ],
    [ 2/3.,   -2/3.,    1/6.,    1/3.,    1/6.,   -2/3.]   ])


    # Coordenadas normalizadas de los nodos 
    # (lo que tenga _ sera indicativo de normalizado, del espacio psi,eta)
    x_, y_ = [None, None], [None, None]
    
    # en [0] el cuadr_8 
    x_[0] = np.array([-1.,  0.,  1., 1., 1., 0., -1., -1.])
    y_[0] = np.array([-1., -1., -1., 0., 1., 1.,  1.,  0.])
    # en [1] el triang_6
    a=r3/2.
    x_[1] = np.array([0.0, 0.5, 1.0, 0.0, -1.0, -0.5])
    y_[1] = np.array([0.0, a,   2*a, 2*a,  2*a,  a])
    
    # Puntos de Gauss y sus pesos
    x_Gaus, y_Gaus, w_Gaus = [],[],[]
    
    # en [0] lo del cuadrado (9 ptos Gauss)
    x_Gaus.append( [-0.774596669241483, 0.000000000000000, 0.774596669241483,-0.774596669241483, 0.000000000000000, 0.774596669241483, -0.774596669241483, 0.000000000000000, 0.774596669241483] )
    y_Gaus.append( [-0.774596669241483,-0.774596669241483,-0.774596669241483, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.774596669241483, 0.774596669241483, 0.774596669241483] )
    w_Gaus.append ( [ 0.308641975308642, 0.493827160493827, 0.308641975308642, 0.493827160493827, 0.790123456790123, 0.493827160493827, 0.308641975308642, 0.493827160493827, 0.308641975308642])

    # en [1] lo del triangulo (6 ptos Gauss)
    x_Gaus.append( [0.000000000000000, 0.000000000000000, 0.725271359470688,-0.337845472747894,-0.725271359470688, 0.337845472747894] )
    y_Gaus.append( [0.317229309127397, 1.544810887650238, 1.573436153005179, 0.959645363743758, 1.573436153005179, 0.959645363743758] )
    w_Gaus.append ([0.190442006391806, 0.386908262797819, 0.190442006391806, 0.386908262797819, 0.190442006391806, 0.386908262797819] )


    # Las triangulaciones para dibujar tensiones etc tambien seran atributos de clase.
    # Cubriran el elemento salvo un margen de 0.01.
    # Aqui tendran los puntos en coordenadas normalizadas (claro). En el programa 
    # hay que llamar a e.dimexy() para poner las reales, y hay  que construir la 
    # malla como adicion de los mallados elementales.
    
    triangulos_=[None, None]
    
    # en [0] lo del cuadrado
    triangulos_[0]=matplotlib.tri.Triangulation(
        [-0.999, 0.999, 0.999, -0.999], [-0.999, -0.999, 0.999, 0.999],
        triangles=[[0,1,2], [2,3,0]])
    triangulos_[0]=matplotlib.tri.UniformTriRefiner(triangulos_[0]).refine_triangulation(subdiv=5)
    
    # en [1] lo del triangulo
    triangulos_[1]=matplotlib.tri.Triangulation(
        [0., 0.999, -0.999], [0.001, 1.73105, 1.73105], 
        triangles=[[0,1,2]])
    triangulos_[1]=matplotlib.tri.UniformTriRefiner(triangulos_[1]).refine_triangulation(subdiv=5)



    def __init__(self, tipoe=-1, nodse=[], xi=[], yi=[]) :
        self.tipoe = tipoe    # tipo de elemento: 0=cuadr_8, 1=triang_6
        self.nodse = nodse    # lista[] con los nodos del elemento
        self.xi = xi
        self.yi = yi    # coordenadas xy de los nodos del elemento


    def dimexy(self, psi,eta):
        # devuelve x,y correspondiente a psi,eta
            # primero calculamos las Hi (psi,eta):
        try:
            unos= np.ones(len(psi)) # por si psi,eta son arrays
        except:
            unos= 1.0
        if   self.tipoe ==0:
            pol=np.array([unos, psi, eta, psi**2, psi*eta, eta**2, psi**2*eta, psi*eta**2])
        elif self.tipoe ==1:
            pol=np.array([unos, psi, eta, psi**2, psi*eta, eta**2])
        Hi= np.matmul(pol.T,self.coefH[self.tipoe])
        
            # ahora las x,y:
        x = np.matmul( Hi, self.xi)
        y = np.matmul( Hi, self.yi )
        return(x,y)
    
    def dimeJaco(self, psi,eta):
        if   self.tipoe ==0:
            pol_psi=np.array([0., 1., 0., 2*psi, eta,   0.,  2*psi*eta,  eta**2])
            pol_eta=np.array([0., 0., 1.,  0.,   psi, 2*eta,   psi**2,  2*psi*eta])
        elif self.tipoe ==1:
            pol_psi=np.array([0., 1., 0., 2*psi, eta,  0.])
            pol_eta=np.array([0., 0., 1.,   0.,  psi, 2*eta])
        Hi_psi= np.matmul(pol_psi, self.coefH[self.tipoe])
        Hi_eta= np.matmul(pol_eta, self.coefH[self.tipoe])
        x_psi= np.matmul(self.xi, Hi_psi)
        x_eta= np.matmul(self.xi, Hi_eta)
        y_psi= np.matmul(self.yi, Hi_psi)
        y_eta= np.matmul(self.yi, Hi_eta)
        jaco=np.array([[x_psi, y_psi],[x_eta, y_eta]]) # mejor doy el determinante
        jacodet= x_psi*y_eta - x_eta*y_psi
        
        return(jacodet)
    
    def pinta_elem(self, colorin='', tipolin='', ancholin=0.8, desplacin=False, uscal=1.):
        # dibuja las lineas (curvas o no) de los lados, nada mas
        npl= 18 # numero de puntos en cada lado, para el trazado
        def puntea_lado(self, xi_, yi_, x, y):
            for i in range(npl):
                a,b= self.dimexy(xi_[i], yi_[i])
                if desplacin:
                    u= np.matmul(N_xy(a,b),a_despl)
                    a += u[0] * uscal
                    b += u[1] * uscal
                x.append(a)
                y.append(b)
            return()

        x, y = [], []
        if self.tipoe == 0:
            xi_ , yi_ = np.linspace(-1,1,npl), [-1]*npl
            puntea_lado(self, xi_, yi_, x, y)
            
            xi_ , yi_ = [1]*npl, np.linspace(-1,1,npl)
            puntea_lado(self, xi_, yi_, x, y)
            
            xi_ , yi_ = np.linspace(1,-1,npl), [1]*npl
            puntea_lado(self, xi_, yi_, x, y)
            
            xi_ , yi_ = [-1]*npl, np.linspace(1,-1,npl)
            puntea_lado(self, xi_, yi_, x, y)
            
        elif self.tipoe==1:
            r3=np.sqrt(3.)

            xi_ , yi_ = np.linspace(-1,1,npl), [r3]*npl
            puntea_lado(self, xi_, yi_, x, y)

            a= np.array([-1, -r3])
            b=np.linspace(0,1,npl)
            xi_, yi_ = [], []
            for i in range(npl):
                v = np.array([1,r3]) + a*b[i]
                xi_.append(v[0])
                yi_.append(v[1])
            puntea_lado(self, xi_, yi_, x, y)

            a= np.array([-1, r3])
            xi_, yi_ = [], []
            for i in range(npl):
                v = a*b[i]
                xi_.append(v[0])
                yi_.append(v[1])
            puntea_lado(self, xi_, yi_, x, y)
        
        a=['-','--'] if tipolin == ''  else [tipolin,tipolin]
        plt.plot(x,y, a[self.tipoe], linewidth=ancholin, color=colorin)
        return()

class ToolTip(object):
    # para los cuadradillos informativos emergentes
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = Label(tw, text=self.text, justify=LEFT,
                      background="#ffffe0", relief=SOLID, borderwidth=1,
                      font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)



# GUI -  primera  ventana (de presentacion y elecciones) 

def presenta_elige():
    global elige,n_f, nfcompleto


    ##### Ventana hija para licencia #####

    def licencia():
        v1=Toplevel(v0)

        texto='''Este programa es Software Libre (Free Software). Como tal se le aplican los términos de la "GNU General Public License", en su versión 2 o bien (como usted prefiera) en una versión posterior. Básicamente, usted puede:

* Usar libremente el programa 
* Realizar copias del mismo y distribuirlas libremente 
* Estudiar su código para aprender cómo funciona 
* Realizar modificaciones/mejoras del programa 

Bajo las siguientes condiciones:  

* Las modificaciones realizadas al programa deben hacerse públicas bajo una licencia como la presente (así el software puede mejorar con aportaciones realizadas 
sobre el trabajo de otros, como se hace en la ciencia).
* Las modificiones y trabajos derivados en general deben incluir el reconocimiento al autor original (no puede decir que es usted quien ha escrito este programa).
* En este caso, debe mencionar al autor original como: Juan Carlos del Caño, profesor en la Escuela de Ingenierías Industriales de la Universidad de Valladolid (Spain)  

Este programa se distribuye con la esperanza de que sea útil, pero SIN NINGUNA GARANTIA, ni siquiera la garantía de comerciabilidad o de adecuación para un propósito 
particular. Lea la GNU General Public License para más detalles.
Usted debería haber recibido una copia de la GNU General Public License junto con este programa. Si no es así, escriba a la Free Software Foundation: 
Inc., 51 Franklin Street, Fifth Floor, 
Boston, MA 02110-1301, USA.'''

        ttk.Label(v1, text='Resumen de términos de licencia', 
            font=('', 16)).grid(row=0, column=0, columnspan=2, pady=5)

        tcaja = Text(v1, width=45, height=30,wrap='word', font=('Sans',9),
            background='#D9D9D9', foreground='green', border=None, padx=20, pady=12)
        tcaja.grid(column=0, row=1, padx=8, sticky=(N,W,E,S))
        tcaja.insert('1.0',texto)

        scb = ttk.Scrollbar(v1, orient=VERTICAL, command=tcaja.yview)
        scb.grid(column=1, row=1, sticky='ns')
        
        tcaja['yscrollcommand'] = scb.set

        ttk.Button(v1, text='Entendido', width=9, command=v1.destroy).grid(
            column=0, row=2, pady=4, columnspan=2)

        tcaja['state']='disabled'

        v1.grid_columnconfigure(0, weight=1)
        v1.grid_rowconfigure(0, weight=1)
        v1.grid_rowconfigure(1, weight=4)
        v1.grid_rowconfigure(2, weight=1)
        v1.geometry('+240+60')

        v1.focus()
        v1.mainloop()


    ##### Ventana hija para breviario #####

    def breviario():

        v2=Toplevel(v0)
            
        texto='''- Para qué sirve este programa:

Realiza una aproximación de Galerkin a un problema elástico bidimensional (tensión o deformación plana). La geometría y las funciones de aproximación son especificadas por el usuario.
El programa calcula el campo de tensiones (sus tres compontes en el plano) y el campo de desplazamientos, y dibuja una gráfica para cada uno de ellos como resultado. 
La salida de texto ofrece valores numéricos más precisos e información adicional.

Se trata de un software de propósito educacional y para investigación, ya que la aproximación en sí de un problema elástico puede en general realizarse con menor esfuerzo 
humano mediante métodos más sistemáticos como el Método de los Elementos Finitos (el cual es por otra parte un caso particular de aproximación de Galerkin).


- Cómo lo hace:

La geometría se especifica mediante elementos rectangulares curvos de ocho nodos y/o elementos triangulares curvos de seis nodos. Los nodos de un elemento deben darse en orden antihorario. 

Es aconsejable que los elementos tengan un área similar, y usar el mínimo número de elementos que sea suficiente para aproximar la geometría. Téngase en cuenta que la aproximación de Galerkin no usa funciones definidas por elementos, y que por otra parte el programa se ocupa internamente de que las integraciones se realicen con suficiente precisión.

La versión actual del programa solamente admite condiciones de contorno homogéneas en desplazamientos (es decir, condiciones de desplazamiento nulo), y no realiza ninguna comprobación acerca de las mismas. Por ello las funciones de aproximación deben elegirse de forma que satisfagan todas ellas dichas condiciones de contorno en desplazamientos. 

Las condiciones de contorno en tensiones pueden ser fuerzas concentradas o distribuídas con forma arbitraria. También se admiten fuerzas de volumen y cargas térmicas, ambas expresadas como funciones en el dominio. Para más detalles se aconseja ejecutar el ejemplo suministrado, al que acompaña un pequeño tutorial.


- Estado de desarrollo del programa:

Esta versión 0.7 del programa es la primera que se hace pública. Las versiones anteriores fueron de desarrollo. Por ello es probable que exista algún 'bug' inadvertido. El autor agradecerá ser informado si se aprecia alguno. Esta versión no implementa aún la interpolación de las condiciones de contorno en desplazamientos, como se ha expuesto. 

Espero que GalerkinAppr le sea útil.
            ____________________________________
'''
        
        ttk.Label(v2, text='Notas breves',
            font=('', 16)).grid(row=0, column=0, columnspan=2, pady=5)

        tcaja = Text(v2, width=45, height=30,wrap='word', font=('Sans',9),
            background='#D9D9D9', foreground='green', border=None, padx=20, pady=12)
        tcaja.grid(column=0, row=1, padx=8, sticky=(N,W,E,S))
        tcaja.insert('1.0',texto)

        scb = ttk.Scrollbar(v2, orient=VERTICAL, command=tcaja.yview)
        scb.grid(column=1, row=1, sticky='ns')
        
        tcaja['yscrollcommand'] = scb.set

        ttk.Button(v2, text='Entendido', width=9, command=v2.destroy).grid(
            column=0, row=2, pady=4, columnspan=2)

        tcaja['state']='disabled'

        v2.grid_columnconfigure(0, weight=1)
        v2.grid_rowconfigure(0, weight=1)
        v2.grid_rowconfigure(1, weight=4)
        v2.grid_rowconfigure(2, weight=1)
        v2.geometry('+250+70')

        v2.focus()
        v2.mainloop()



    ##### Para el prb por defecto #####

    def default_prb():
        global elige
        elige='default' # sin mas ventanas
        v0.destroy()



    ##### La ventana de inicio (por fin) #####

    v0.title("GalerkinAppr v0.7")

    cuadro = ttk.Frame(v0, padding='9 3 3 3') 
    cuadro.grid(column=0, row=0, sticky=(N, W, E, S))

    ttk.Label(cuadro, text='GalerkinAppr v0.7', font=('', 40)).grid(row=0,
        column=0, columnspan=4)
    ttk.Label(cuadro, text='Aproximador de Galerkin para problemas 2D', 
        font=('Courier', 16)).grid(row=1, column=0, columnspan=4)
    ttk.Label(cuadro, text='by:   Juan Carlos del Caño\n').grid(row=2,
        column=0, columnspan=4)

    # hago la parte izda de la ventana

    ttk.Separator(cuadro, orient=HORIZONTAL).grid(column=0, row=3,
        columnspan=4, sticky='ew')
        
    texto= 'Esto es Software Libre (Free Software), lo cual\n'
    texto +='le otorga a usted algunos derechos pero también\n'
    texto +='conlleva algunas obligaciones. \nPor favor lea la licencia.'
    ttk.Label(cuadro, text=texto, foreground='green').grid(row=4, 
        column=0, columnspan=3, sticky='w')

    ttk.Separator(cuadro, orient=HORIZONTAL).grid(row=5, column=0,
        columnspan=3, sticky='ew')
        
    texto=  'Si usted no conoce para qué sirve este programa,\n'
    texto +='o tiene una duda básica respecto del mismo,\n'
    texto +='consulte estas breves notas.\n'
    texto +='O bien consulte la documentación.'
    ttk.Label(cuadro, text=texto, foreground='green').grid(row=6, 
        column=0, columnspan=3, sticky='w')

    ttk.Separator(cuadro, orient=HORIZONTAL).grid(row=7, column=0,
        columnspan=3, sticky='ew')
        
    texto = 'Puede empezar un problema nuevo desde cero.\n'
    texto +='Deberá especificar la geometría de la sección\n'
    texto +='(puntos y tramos) y algunos otros datos generales\n'
    texto +='del problema. También puede cargar desde archivo\n'
    texto +='un problema que haya sido guardado previamente.'
    ttk.Label(cuadro, text=texto, foreground='green').grid(row=8,
        column=0, columnspan=3, sticky='w')

    ttk.Separator(cuadro,orient=HORIZONTAL).grid(row=9,column=0,sticky='ew')
    ttk.Separator(cuadro,orient=HORIZONTAL).grid(row=9,column=2,sticky='ew')
    ttk.Label(cuadro, text='o bien:').grid(row=9,column=1)
        
    texto = 'Puede cargar un ejemplo preprogramado internamente\n'
    texto +='que ilustra las capacidades básicas del programa,\n'
    texto +='y con el que puede comprobarse si la instalación de\n'
    texto +='GalerkinAppr fue correcta.'

    ttk.Label(cuadro, text=texto, foreground='green').grid(row=10,
        column=0, columnspan=3, sticky='w')
        
    ttk.Separator(cuadro, orient=HORIZONTAL).grid(row=11, column=0,
        columnspan=4, sticky='ew')


    # ahora hago la parte derecha

    ttk.Separator(cuadro,orient=VERTICAL).grid(row=3,column=3,
        rowspan=9, sticky='ns')
    ttk.Button(cuadro, text='Licencia',
        command=licencia).grid(row=4, column=3)

    ttk.Button(cuadro, text='Breviario', command=breviario).grid(
        row=6, column=3)

    ttk.Button(cuadro, text='Resolver', command=v0.destroy).grid(
        row=8, column=3)

    ttk.Button(cuadro, text='Ejemplo', command=default_prb).grid(
        row=10, column=3)

    for hijo in cuadro.winfo_children():
        hijo.grid_configure(padx=12, pady=8)

    v0.geometry('+70-70')
    v0.focus()
    v0.mainloop()
    
    return()

###################### Fin de ventana de inicio #####################

def rellena_default():
        # inspirado en galerkin_04.dat. Se trata de rellenar como si se
        # hiciese a mano. Aqui todo es texto por tanto.
    global filas_nodos, filas_elems, entry_young, entry_poiss, entry_dilat
    global filas_funcs, filas_fuerzas, filas_xraya, filas_xgrande, filas_temperatura 
    global nnodos, nelems, nff, tens_plana  # este es StringVar()
    
    nnodos=22
    lista= [
    '1',    ' 0.000',   '-2.500',
    '2',    ' 2.500',   '-2.500',
    '3',    ' 5.000',   '-2.500',
    '4',    '10.000',   '-2.500',
    '5',    '15.000',   '-2.500',
    '6',    '17.500',   '-2.500',
    '7',    '20.000',   '-2.500',
    '8',    ' 0.000',   '0.0',
    '9',    ' 1.900',   '0.0',
   '10',    ' 8.200',   '0.0',
   '11',    '12.000',   '0.0',
   '12',    '14.000',   '0.0',
   '13',    '20.000',   '0.0',
   '14',    ' 0.000',   '2.500',
   '15',    ' 2.400',   '2.500',
   '16',    ' 5.000',   '2.500',
   '17',    ' 7.600',   '2.500',
   '18',    '10.000',   '2.500',
   '19',    '12.400',   '2.500',
   '20',    '15.000',   '2.500',
   '21',    '17.600',   '2.500',
   '22',    '20.000',   '2.500' ]
    
    for i in range(nnodos):
        filas_nodos[i][0].insert(0,lista[3*i])
        filas_nodos[i][1].insert(0,lista[3*i+1])
        filas_nodos[i][2].insert(0,lista[3*i+2])

    nelems=5
        # ielem, nodos  (el tipoe es una cosa interna)
    lista=[
    ['1' , '1', '9' ,'16' , '15', '14', '8'] ,
    ['2' , '1', '2' ,'3'  , '10', '18', '17', '16', '9'],
    ['3' , '3', '4' ,'5'  , '11', '18', '10'],
    ['4' , '5', '12' ,'20', '19', '18', '11'],
    ['5' , '5', '6' ,'7'  , '13', '22', '21', '20', '12']
    ]
    
    for i in range(nelems):
        for j in range (len(lista[i])):
            filas_elems[i][j].insert(0,lista[i][j])
    
    
    # ctes elasticas
    entry_young.insert(0,'1000.')
    entry_poiss.insert(0,'0.')
    entry_dilat.insert(0,'0.')
    tens_plana.set('1')



    # funciones de aproximacion
    nff=9
    lista= ['x* y       /1.e2',
    'x* y**3    /1.e4',
    'x* y**5    /1.e6',
    'x**2 *y    /1.e3',
    'x**2* y**3 /1.e5',
    'x**2* y**5 /1.e7',
    'x**3* sin(y)/1.e3',
    'x**3* y**3 /1.e6',
    'x**3* y**5 /1.e8' ]  # eso es para ux
    for i in range(nff):
        filas_funcs[i][1].insert(0,lista[i].strip())

    lista=['1. -exp(-x)',
    'x    *y**2  /1.e3',
    'x    *y**4  /1.e5',
    'x**2        /1.e2',
    'x**2 *y**2  /1.e4',
    'x**2 *y**4  /1.e6',
    'x**3        /1.e3',
    'x**3 *y**2  /1.e5',
    'x**3 *y**4  /1.e7' ] # eso es para uy
    for i in range(nff):
        filas_funcs[i][2].insert(0,lista[i].strip())

    # cargas
    '''
    filas_xraya[0][0].set('1')
    filas_xraya[0][1].insert(0,'17')
    filas_xraya[0][2].insert(0,'18')
    filas_xraya[0][3].insert(0,'19')
    filas_xraya[0][4].insert(0,'0.')
    filas_xraya[0][5].insert(0,'-0.3125') # es ~aprox a F=-1 en n18 (osea x=10)
    filas_xraya[0][6].insert(0,'0.')
    filas_xraya[0][7].insert(0,'0.')
    filas_xraya[0][8].insert(0,'0.')
    filas_xraya[0][9].insert(0,'0.')
    '''
    
    filas_fuerzas[0][0].insert(0,'22')
    filas_fuerzas[0][1].insert(0,'0.0')
    filas_fuerzas[0][2].insert(0,'-0.25') # en el gui son columnas, en [] son filas

    filas_fuerzas[1][0].insert(0,'13')
    filas_fuerzas[1][1].insert(0,'0.0')
    filas_fuerzas[1][2].insert(0,'-0.5')

    filas_fuerzas[2][0].insert(0,'7')
    filas_fuerzas[2][1].insert(0,'0.0')
    filas_fuerzas[2][2].insert(0,'-0.25')
    
    '''
        # lo siguiente es para probar, no estara en ej final
    filas_xgrande[0].insert(0,'x*(x-10.)*(y-4)')
    filas_xgrande[1].insert(0,'x*y*y')
    filas_temperatura[0].insert(0,'x*(x-2.)*(x-4.)')
    '''

def rollete_default():
    v7=Toplevel(v0)
    
    texto= '''El ejemplo considera un sólido de geometría rectangular de 5 x 20 en el sistema de unidades que el usuario haya elegido y que el programa asume de forma transparente. Los cálculos se realizan para una profundidad igual a uno en esas mismas unidades.


- Pestaña 'Dominio' -

- En esta pestaña se especifican la geometría y datos físicos del sólido. En este caso la geometría se ha discretizado  en 5 elementos de lados parabólicos derivados bien sea de un cuadrado de 8 nodos o bien de un triángulo equilátero de 6 nodos. Es aconsejable que los ángulos de los elementos no estén excesivamente distorsionados respecto de su geometría canónica.

- Las coordenadas de los nodos están referidas a un eje 'x' horizontal situado a mitad de la altura del rectángulo, y a un eje 'y' vertical situado sobre el lado izquierdo. La numeración de nodos no tiene porqué ser consecutiva ni completa, aunque se suela hacer así por claridad cuando no hay un motivo en contra.

- Los elementos se especifican mediante un número identificativo (nuevamente esta numeración no tiene porqué ser consecutiva ni completa), seguido de sus nodos numerados consecutivamente en sentido antihorario. La numeración debe comenzar en un vértice del elemento. Pulsando el botón 'Visualizar' del contexto 'Elementos' se muestra una representación de los datos introducidos. Ello es de interés para comprobación de datos cuando los mismos han sido introducidos a mano. 

- El programa identifica el tipo de elemento a partir del número de nodos que el usuario haya introducido (seis para triangular, ocho para rectangular). Por ello es importante no dejar caracteres espúreos en casillas que debieran estar vacías. Lo mismo se aplica a las funciones de aproximación y a las cargas de contorno y de volumen de los apartados siguientes.

- Recuérdese que la aproximación de Galerkin NO usa funciones definidas sobre 'elementos' en los cálculos, siendo de hecho extraño a este método el propio concepto de elemento. Independientemente de ello, es necesario describir de alguna manera el dominio de integración (el ocupado por el sólido). Eventualmente, este programa utiliza para ello un enfoque basado en elementos aunque bien podría haberse adoptado algún otro enfoque. Con lo anterior se intenta resaltar que no se obtendrá una aproximación mejor por usar más elementos, y que tiene sentido usar solamente el menor número de elementos que sea suficiente para aproximar razonablemente la geometría (y en su caso para que haya nodos suficientes para especificar las cargas de contorno). El ejemplo tiene fines ilustrativos y por ello se ha complicado ese aspecto usando más elementos de los necesarios e incluso curvando sus lados. En un uso típico hubiese bastado con un solo elemento rectangular para todo el sólido.

- Finalmente en la pestaña 'Dominio' se pueden especificar las constantes elásticas del material y el tipo de problema bidimensional (tensión o deformación plana). En el ejemplo figura a título ilustrativo un módulo de Young E=1000, y coeficientes de Poisson y de dilatación nulos, independientemente de lo realistas que estos valores pudieran ser. 


- Pestaña 'Funciones' -

- Este es el lugar donde especificamos las funciones de aproximación cuya combinación construirá el campo de desplazamientos aproximado. Debe haber el mismo número de funciones para aproximar las componentes 'ux' y 'uy' de desplazamiento. 

- Por motivos de seguridad las entradas de funciones solo admiten 'x' e 'y' como variables, los dígitos 0 a 9, los cuatro signos aritméticos, el punto decimal, la 'e' para exponentes de 10 en constantes, los paréntesis y las palabras especiales 'exp' (exponenciación de base el número e), 'sin' (seno de un ángulo en radianes), 'cos' (ídem. para el coseno) y 'log' (logaritmo natural). La exponenciación general se denota con **. Obsérvese que la mayoría de estas funcionalidades se han utilizado en el ejemplo con fines ilustrativos.

- En su versión actual el programa sólo admite condiciones de contorno homogéneas en desplazamientos (del tipo desplazamiento nulo en una o en sus dos componentes), y las mismas deben imponerse a través de las funciones de aproximación, las cuales deben anularse en los puntos o líneas correspondientes. En el ejemplo los puntos del lado izquierdo tienen impedido el movimiento en ambas direcciones, por lo que tanto las funciones de 'ux' como las de 'uy' tienen algún factor que se anula en x=0 ('x' elevado a alguna potencia u otra función de 'x'). 

- Apréciese que todas las funciones han sido "normalizadas" dividiendo por una longitud representativa del problema (se ha tomado 10 en este caso) elevado a la potencia del término polinómico en xy presente en la función. Esto no es necesario pero facilita la comparación directa de los coeficientes de la aproximación por obtener, para juzgar qué funciones han tenido más peso (fueron más "acertadas") etc. 

- El "acertar" con una buena aproximación al campo de desplazamientos usando pocas funciones (una o dos) es bastante difícil salvo para casos triviales. Usando un número moderado de funciones (digamos 9 como en el ejemplo) el pronóstico es mucho mejor, especialmente si las funciones han sido elegidas con algún cuidado. En nuestro caso las cargas consistirán básicamente en una fuerza vertical en el extremo derecho del sólido, lo que implica que el problema será antisimétrico. Las funciones de aproximación se han elegido impares en 'y' para ux, y pares en 'y' para uy, lo que se adapta a la antisimetría mejorando con ello las expectativas de la aproximación.


- Pestaña 'Cargas' -

- En el ejemplo se consideran tres cargas puntuales hacia abajo en los nodos 7, 13 y 22, de valores 0.25 , 0.5 , 0.25 respectivamente, en el sistema de unidades coherente que el usuario estuviese empleando. Apréciese cómo están especificadas en la interfaz. Los otros tipos de carga no actúan en nuestro ejemplo pero se comentarán brevemente más adelante. Deben quedar en blanco por ahora.


- Análisis de resultados -

- La geometría, sustentación y cargas de nuestro sólido recuerdan a una viga corta empotrada y con una carga unidad en la punta. La solución clásica de Resistencia de Materiales (modelo de Navier-Bernuilli, aplicado "cerrando los ojos") predice un desplazamiento en la punta de 0.256 unidades de desplazamiento, y unas tensiones máximas en el empotramiento (nodos 1 y 14 en nuestro caso) de 4.8 unidades de fuerza por unidad de área. Usaremos estos resultados como referencia comparativa.

- La salida de texto indica un desplazamiento uy prácticamente igual en los nodos 7, 13 y 22, de valor 0.265. El resultado compara bien con el de la Resistencia de Materiales, y es de esperar que sea más exacto (¿porqué?). Las tensiones en los nodos 1 y 14 figuran como 4.74 (salvo el signo) lo que nuevamente compara bien con el modelo de Resistencia de Materiales. 

- Como observación previa acerca de las ventanas gráficas: no respetan la proporción ancho/alto del sólido para aprovechar mejor el área de trazado, pero si ud. modifica manualmente la geometría de la ventana la figura le acompañará. Existe una opción con la que puede imponer que se respeten las proporciones reales, pero mejor no distraer la atención ahora con este detalle.

- La gráfica de sigma_xx muestra que esta componente de tensión es pequeña en el extremo derecho y toma valores máximos mayores según nos movemos hacia la izquierda, de tracción por arriba y de compresión por abajo. En el eje x se mantiene nula. Todo ello es grosso modo lo esperado para el funcionamiento de una viga en cuanto a esta componente de tensión, que es típicamente la más relevante. 

- La gráfica de sigma_yy muestra valores muy próximos a cero en todo el dominio salvo en las esquinas con carga aplicada, que nuevamente es grosso modo lo que cabe esperar. 

- Esta componente sigma_yy es en realidad nula en las superficies superior e inferior del sólido (omitiendo las esquinas), pero no se obtiene así porque los métodos numéricos aproximados, como el Método de los Elementos Finitos y la propia Aproximación de Galerkin, aproximan estas condiciones de contorno de forma compatible con las funciones de aproximación utilizadas, no de forma exacta. Como observación final véase que en el eje x se obtiene sigma_yy exactamente nulo, exactitud que deriva de que hemos usado funciones de desplazamiento exactamente antisimétricas.

- La gráfica de sigma_xy muestra como tendencia general una concentración mayor junto al eje 'x', mientras que las líneas de sigma_xy=0 se mantienen próximas a los contornos superior e inferior. El valor de estas tensiones es de casi dos órdenes de magnitud menor que sigma_xx, y son negativas como cabe esperar. El modelo de la Resistencia de Materiales predice de hecho una distribución parabólica en 'y', con el máximo en y=0 y que se anula en los contornos superior e inferior. Esa distribución parabólica sería constante en 'x' según dicho modelo. 

- Para que lo anterior quede cuantificado, observemos en la salida de texto los valores de sigma_xy para los nodos 8, 9, 10, 11, 12 y 13, en los que esperamos el máximo. Apreciamos que dichos valores están en un estrecho intervalo en torno a 0.31, salvo el nodo 9. El valor que predice el modelo de Resistencia de materiales es de 0.30.

- La gráfica de desplazamientos muestra la configuración deformada (con desplazamientos exagerados) superpuesta a la configuración inicial, y un diagrama de flechas orientativo para el campo de desplazamientos. Se trata de una figura que permite apreciar a grandes rasgos la naturaleza de la solución obtenida, e identificar posibles errores en los datos etc. 

- Con el fin de usar esta figura para algo más, podemos ampliar la deformada del contorno derecho (la herramienta lupa en la parte inferior de la ventana permite hacerlo) para apreciar el alabeo de la sección. No apreciamos curvatura, por lo que concluimos que esta sección debe alabear muy poco. Esto es coherente con nuestros conocimientos de Resistencia de Materiales según los cuales las secciones macizas alabean muy poco. También podríamos pensar que el modelo que hemos propuesto no ha sido capaz de recoger dicho efecto, pero esto es poco probable ya que hemos usado un buen número de funciones de desplazamiento 'ux' con potencias impares de 'y' diversas. 


- Variaciones propuestas -

- Puede pulsar el botón "cerrar gráficos" de la interfaz principal en este momento. Es preciso hacerlo cuando se van a modificar datos y recalcular el problema, como es el caso. También cuando se va a resolver un problema nuevo.

- Volviendo a la pestaña 'Cargas', borremos las casillas de cargas puntuales, cuidando de no dejar espacios en blanco espúreos etc. Podemos hacerlo navegando por cada casilla ocupada (las teclas de tabulación y de borrar permiten hacerlo rápidamente) o bien utilizar el botón "limpiar" y eligiendo "cargas". Añadimos después una carga distribuida de contorno que abarcará los nodos 17, 18 y 19 , siendo nula en el primero y el último y máxima en el 18. Este valor máximo se tomará como 0.3125  hacia abajo para que la resultante de la distribución sea igual a uno. En la primera fila de "Cargas distribuidas de contorno" marcamos el pulsador de la izquierda, con lo que estamos indicando que especificaremos las componentes de los vectores en coordenadas normal y tangencial ('nt') en lugar de hacerlo en las coordenadas globales 'xy'. El resto de la fila se rellena con los siguientes datos:

  19  18  17     0.  -0.3125  0.       0.  0.  0.

- La convención es: primero se escriben los tres nodos en en orden de avance elegido sobre le contorno, después las componentes 'x' (o bien normales) en cada uno de esos nodos, y después las componentes 'y' (o bien tangenciales). La normal está dirigida hacia la derecha según avance de los nodos, mientras que la tangencial sigue el sentido de avance de los nodos, en el orden en que se han proporcionado. El orden 19-18-17 se ha adoptado para que la normal sea la exterior al sólido (como se hace para los elementos) pero no daría error el adoptar 17-18-19 y cambiar de signo las componentes de vector tensión suministradas. 

- Nótese también que la región de carga abarca parcialmente dos elementos. Esto está permitido porque el programa no usa internamente los elementos para tratar este tipo de cargas. Por lo demás, la implementación realizada describe estas cargas de contorno por tramos parabólicos como el mostrado. Deben preverse nodos adecuados en la discretización de la geometría para poder especificar las cargas de este tipo que haya en el problema. Evidentemente, un distribución de carga más complicada deberá describirse con varios tramos parabólicos. Una carga lineal ó constante es para el programa un caso particular de lo anterior y requiere igualmente tres nodos del contorno para ser descrita.

- Pulsando el botón "calcular" se rehacen los cálculos con los nuevos datos. Se aprecia que la evolución de las tensiones sigma_xx tiene la forma esperada para una viga en la mitad izquierda, mientras que en la mitad derecha son significativamente pequeñas (como es de esperar ya que sería una "zona descargada" de la viga). En la gráfica de sigma_xy la linea de valor nulo "intenta" adaptarse al contorno, salvo evidentemente al izquierdo, y las tensiones sigma_yy se mantienen próximas a cero en la parte del dominio en la que se espera que así sea.

- Merece la pena comentar que en la gráfica de sigma_yy se obtienen tensiones algo mayores en la zona de carga. Esto era de esperar ya que la propia carga es localmente una tensión sigma_yy que el modelo intentará ajustar. Lo que puede sorprender a primera vista es que en la parte inferior (descargada) aparezca una evolución similar y de signo contrario. Esto se debe a que habíamos ajustado las funciones de aproximación para describir un problema antisimétrico, y éste ya no lo es. Si representásemos las tensiones obtenidas en el contorno encontraríamos algo parecido a la mitad de la carga asignada al contorno superior (y de compresión) y la otra mitad al contorno inferior (y de tracción), formando así un problema genuinamente antisimétrico de alguna manera parecido al original.

- Para justificar que el resultado es razonable podemos hacer una comparación "quick'n'dirty" con una viga de longitud 20 con carga unidad en su mitad. El desplazamiento en la punta y las tensiones sigma_xx máximas en el empotramiento son:
            Modelo de RM    GalerkinAppr
u_máx          0.080           0.084
tens_máx       2.4             2.0
Que teniendo en cuenta las varias inexactitudes asumidas a sabiendas (viga ya ultra-corta de relación 2:1, principio de Saint-Venant mal aplicado, problema no antisimétrico), parece un resultado razonable al menos en cuanto a orden de magnitud.

- Exploremos ahora las cargas de volumen. Se expresan mediante dos funciones de las coordenadas espaciales, las cuales tendrán dimensiones de fuerza por unidad de volumen: una función para la componente 'x' de la fuerza de volumen y otra para la componente 'y'. Para mantener la antisimetría probemos con una carga de dirección 'y' dada por la función -0.01*x (el 0.01 tendría dimensiones: fuerza partido por longitud a la 4). Borramos cualquier otra carga, cerramos las ventanas gráficas y ponemos 0.0 en la componente 'x' & -0.01*x en la componente 'y' de las cargas de volumen. Pulsamos el botón "calcular".

- Nuevamente cabe aproximar la configuración por una viga (a sabiendas de que es corta)  según el modelo de Resistencia de Materiales. En este caso la función unidimensional de cargas p_y de dicho modelo será 0.01*x*5*1 = 0.05*x hacia abajo, donde 0.05 tiene dimensiones de fuerza partido por longitud al cuadrado. La comparación de este modelo con los resultados obtenidos es:
            Modelo de RM    GalerkinAppr
u_máx          1.41            1.46
tens_máx       32.0            30.0
Comparación que nuevamente anima a pensar que los cálculos parecen correctos dentro de la aproximación planteada.

- Recalculemos para una carga de volumen simétrica, digamos que 0.01*x en la componente 'x', cero en la componente 'y'. Cerramos las gráficas, ponemos estos valores en el apartado de cargas de volumen, y pulsamos el botón fx<>fy para intercambiar las funciones, que ahora corresponderán a un problema simétrico. Podemos modificar también las funciones de aproximación, lo que haremos en este momento cambiando la exponencial de la primera función de ux (que se puso simplemente como muestra) por la sencilla x/10. El botón calcular genera los resultados. Podemos comparar con los valores teóricos, que se obtienen fácilmente debido a que el coeficiente de Poisson es nulo:
            Teórico      GalerkinAppr
u_máx       0.02666        0.02666
tens_máx    2.00           1.9999 
Se aprecia que en este caso sencillo las funciones elegidas son sensiblemente capaces de reproducir la solución exacta. Las gráficas muestran la evolución esperada de sigma_xx (constante en 'y', máxima en x=0), mientras que las otras componentes de tensión tienen valores muy próximos al cero teórico.
Si lo desea, recalcule para un coeficiente de Poisson 0.3 y aprecie las diferencias en las gráficas, especialmente en el entorno de x=0.

- Finalmente ensayemos un campo de temperaturas. Sea 25 la temperatura a la izquierda y 5 a la derecha variando linealmente. Esto es T= (25-x). Nótese que el problema es nuevamente simétrico por lo que utilizaremos las mismas funciones que en el caso anterior. Cerramos por tanto las ventanas gráficas abiertas, borramos cualquier otra carga, ponemos (25-x) como carga térmica, en la pestaña 'Dominio' ponemos un coeficiente de dilatación de 1.0e-3, y finalmente pulsamos el botón "calcular". 

- Podemos estimar el desplazamiento máximo como el incremento de longitud de la viga si no tuviese sustentación, y las tensiones máximas (que serán sigma_yy) con las asociadas a la dilatación vertical impedida a la izquierda. La comparación resulta:
            Esperado      GalerkinAppr
u_máx        0.300          0.303
tens_máx     25.0           25.0 
Con coincidencia incluso mejor de lo que cabía esperar para problemas que sólo 'se parecen'. 

- La gráfica de desplazamientos muestra un engrosamiento local (vertical) que crece aproximadamente desde la cota del nodo 14 hasta la del 17 para luego decrecer. Esto es esperable dado que es una zona de temperaturas aún altas, que van decreciendo hacia la derecha. La gráfica de sigma_yy muestra tensiones bajas en toda la longitud salvo en la zona de la izquierda, como cabía esperar. 

- Para finalizar cabe apuntar que si nuestro problema no es simérico ni antisimétrico, bien sea debido a las cargas bien sea porque la geometría no es simétrica, no podremos apoyarnos en este criterio para elegir funciones de aproximación. De hecho será adecuado incluir funciones que satisfagan tanto la simetría como la antisimetría, así como funciones generales. El criterio será simplemente que parezca probable que entre las combinaciones posibles de esas funciones queden cubiertas las posibilidades de movimiento del sólido que parezcan más relevantes. Y por supuesto debe evitarse incluir funciones linealmente dependientes de otras ya presentes, y demás aspectos relacionados con la manera de operar de la Aproximación de Galerkin. Dichos aspectos exceden el propósito de esta breve descripción de uso del programa, pero es particularmente necesario haber comprendido las bases del método (¡haberlo estudiado!) para realizar un uso correcto del mismo.
    '''



    ttk.Label(v7, text='Aprendiendo del ejemplo.', background='#F0F0F0',
        font=('', 16)).grid(row=0, column=0, columnspan=2, pady=5)

    tcaja = Text(v7, width=50, height=35,wrap='word', font=('Sans',9),
        background='#F0F0F0', foreground='green', border=None, padx=20, pady=12)
    tcaja.grid(column=0, row=1, padx=8, sticky=(N,W,E,S))
    tcaja.insert('1.0',texto)

    scb = ttk.Scrollbar(v7, orient=VERTICAL, command=tcaja.yview)
    scb.grid(column=1, row=1, sticky='ns')
    
    tcaja['yscrollcommand'] = scb.set

    ttk.Button(v7, text='Entendido', width=9, command=v7.destroy).grid(
        column=0, row=2, pady=4, columnspan=2)

    tcaja['state']='disabled'

    v7.grid_columnconfigure(0, weight=1)
    v7.grid_rowconfigure(0, weight=1)
    v7.grid_rowconfigure(1, weight=4)
    v7.grid_rowconfigure(2, weight=1)
    v7.geometry('+640+90')



def obviar_gui():
    global elige, n_f, nfcompleto, nodos, elems, nnodos, nelems
    global ctes, nff, ffx, ffy, cargas
    global extremosx, extremosy, anchox, anchoy, sin_gui
    

    texto ='''Esta función permite resolver problemas con más nodos, elementos, funciones, y cargas de lo que puede acomodar la interfaz. Si continúa se procesará directamente el fichero de datos que indique (y que ud. habrá preparado típicamente a mano). Si todo va bien aparecerán directamente los resultados del análisis.\nLas opciones que haya especificado se respetarán, vuelva atrás si quiere modificarlas.\nCuando acabe cierre los gráficos y podrá seguir usando la interfaz normalmente. 
    '''
    sin_gui=messagebox.askokcancel(title='Por favor confirme', message='¿Sabía que...?',
        detail=texto,icon='question', parent=v0)
    if not sin_gui: return()
    
    a_borrar() # que no se vean cosas anteriores que no son
    a_cargar()
    
    if valida_ff(ffx) or valida_ff(ffy): return()
    if cargas.volumen:
        if valida_ff(cargas.volumen): return()
    if cargas.termica:
        if valida_ff(cargas.termica): return()
    
    extremosx, extremosy = [0.,0.], [0.,0.] # cada uno [min, max]
    for n in nodos.values():
        if n.x < extremosx[0]: extremosx[0]=n.x
        if n.x > extremosx[1]: extremosx[1]=n.x
        if n.y < extremosy[0]: extremosy[0]=n.y
        if n.y > extremosy[1]: extremosy[1]=n.y
    anchox = extremosx[1]-extremosx[0]
    anchoy = extremosy[1]-extremosy[0]
        
    
    ver_ff()
    calcula()
    #sin_gui = False  # libera el gui por si se quieren hacer mas cosas

    return()



def valida_ff(ff):
    # entrada= una lista de funciones de forma, o funciones de carga...
    bien=' 0123456789+-*/.exy()'
    for ffi in ff:
        a=ffi
        a.strip()
        esta_N=a
        a= a.replace('exp','').replace('sin','').replace('cos','').replace('log','')
        for b in bien:
            a=a.replace(b,'')
        if len(a): # no debería quedar nada
            texto= 'No se entiende la expresión: \n'
            texto+= esta_N +'\ncomo una función válida.'
            texto+='\nPor favor revísela.'
            messagebox.showinfo(message=texto,title='MAL',parent=v0)
            return(1)
        else:
            try:
                a=esta_N
                x,y= 2.2, 3.3
                kk=eval(esta_N)
            except (TypeError, ValueError, SyntaxError, NameError, AttributeError):
                texto= 'No se puede evaluar la función:\n'
                texto += esta_N or '   - (casilla en blanco) -'
                texto+='\nRevise su sintaxis.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v0)
                return(1)
    return(0)

def flip_f():
    for i in range (nff):
        ffx, ffy = filas_funcs[i][1].get(), filas_funcs[i][2].get()
        filas_funcs[i][1].delete(0,'end')
        filas_funcs[i][2].delete(0,'end')
        filas_funcs[i][1].insert(0,ffy)
        filas_funcs[i][2].insert(0,ffx)





def niveles_color():
    global niveles_sigma, niveles_sigma_bak # puede venir (y salir) con blancos
    
    def extremosNiveles():
        for i in range(7):
            filas_niveles[i][0].delete(0,'end')
            filas_niveles[i][1].delete(0,'end')
            filas_niveles[i][0].insert(0,str(niveles_sigma_bak[i][0]) )
            filas_niveles[i][1].insert(0,str(niveles_sigma_bak[i][1]))
    
    def redibujarNiveles():
        for i in range(7):
            a,b = filas_niveles[i][0].get() , filas_niveles[i][1].get()
            if a=='' or b=='':
                niveles_sigma[i][0],niveles_sigma[i][1] = '', ''
            else:
                try:
                    niveles_sigma[i][0], niveles_sigma[i][1]= float(a), float(b)
                except (ValueError, SyntaxError):
                    texto ='Hay un error en la fila '+str(i+1)
                    texto+=' Por favor compruébelo.'
                    messagebox.showinfo(parent=v6,message= texto)
                    return()
                if niveles_sigma[i][0] > niveles_sigma[i][1]:
                    texto ='El nivel min debe ser menor que' 
                    texto+='el max. En la fila '+str(i+1)+' esto'
                    texto+='\nno se cumple.'
                    messagebox.showinfo(parent=v6,message= texto)
                    return()
        for i in plt.get_fignums(): plt.close(i)
        salida_grafica()
        v6.focus()

    def copiarNiveles():
        kk1, kk2 = filas_niveles[0][0].get(), filas_niveles[0][1].get()
        for i in range(1,7):
            filas_niveles[i][0].delete(0,'end')
            filas_niveles[i][1].delete(0,'end')
            filas_niveles[i][0].insert(0,kk1)
            filas_niveles[i][1].insert(0,kk2)
    
    v6 =Toplevel(v0)
    v6.title('Niveles de color')
    texto = 'Seguidamente puede especificar el rango de valores\n'
    texto+= 'que será representado por colores en las gráficas\n'
    texto+= 'de tensiones. Por ejemplo si el problema tiene una\n'
    texto+= 'carga puntual, probablemente le interese poco lo\n'
    texto+= 'que el método haga para aproximar la singularidad\n'
    texto+= 'localmente, siendo de más interés lo que ocurre en el\n'
    texto+= 'resto del sólido.\n'
    texto+= 'Si deja alguna casilla en blanco la componente de\n'
    texto+= 'tensión correspondiente se dibujará en todo su rango,\n'
    texto+= 'tal como se obtiene tras "calcular".\n'
    texto+= 'Utilícese como segundo paso tras "calcular" (para\n'
    texto+= 'poder juzgar el rango de interés).\n'
    ttk.Label(v6, text=texto, background='#EDECEB', foreground='green').grid(
            row=0, column=0, columnspan=5, padx=8, pady=6)
    ttk.Label(v6, text='min', background='#EDECEB').grid(row=1,column=1)
    ttk.Label(v6, text='max', background='#EDECEB').grid(row=1,column=2)
    ttk.Label(v6, text='s_xx', background='#EDECEB').grid(row=2,column=0)
    ttk.Label(v6, text='s_yy', background='#EDECEB').grid(row=3,column=0)
    ttk.Label(v6, text='s_xy', background='#EDECEB').grid(row=4,column=0)
    ttk.Label(v6, text='s_vM', background='#EDECEB').grid(row=5,column=0)
    ttk.Label(v6, text='s_rr', background='#EDECEB').grid(row=6,column=0)
    ttk.Label(v6, text='s_tt', background='#EDECEB').grid(row=7,column=0)
    ttk.Label(v6, text='s_rt', background='#EDECEB').grid(row=8,column=0)
    ttk.Separator(v6, orient=HORIZONTAL).grid(row=2, column=3,sticky='ew')
    
    filas_niveles=[]
    for i in range(7):
        e1= ttk.Entry(v6, width=8)
        e1.grid(row=i+2, column=1)
        e1.insert(0,str(niveles_sigma[i][0]))
        
        e2= ttk.Entry(v6, width=8)
        e2.grid(row=i+2, column=2)
        e2.insert(0,str(niveles_sigma[i][1]))
        filas_niveles.append([e1,e2])
    
    ttk.Button(v6, text='obtener rangos', command=extremosNiveles).grid(
        row=9, column=0, columnspan=2, padx=[0,10], sticky='e')
    ttk.Button(v6, text='redibujar', command=redibujarNiveles).grid(
            row=9, column=2, columnspan=3, padx=[10,0], sticky='w')
    ttk.Button(v6, text='c\no\np\ni\na\nr\n \ne\nn\n \nt\no\nd\no\ns', width=2,
            command=copiarNiveles).grid(row=2, column=4, rowspan=7, sticky='w')

    for hijo in v6.winfo_children(): hijo.grid_configure(padx=4,
        pady=3)
    v6.focus()
    v6.mainloop()


def obten_datos():
    global nnodos,nelems,nff, nodos,elems,ffx,ffy, ctes,cargas 
    global extremosx, extremosy, anchox, anchoy
    
    nodos={}
    for fila in filas_nodos:
        n0=Nodo(x=0., y=0., ux=0., uy=0., elemn=[], sigma=[0.,0.,0.])
        x0,x1,x2=fila[0].get(), fila[1].get(), fila[2].get()
        if x0 != '':
            try:
                i=int(x0)
                x=float(x1)
                y=float(x2)
            except (TypeError, ValueError):
                texto= 'Hay un error en los nodos[]. Compruebe'
                texto+=' que están todos los datos, etc.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v0)
                return(1)
        n0.x, n0.y = x,y
        nodos[i]= n0
        nnodos=len(nodos)   # mas seguro actualizar a cada paso

    elems={}
    for fila in filas_elems:
        e0=Elem(tipoe=-1, nodse=[], xi=[], yi=[]) # los tipos mutables requieren
                                                  # inicializacion explicita
        lista_txt=[]
        for i in range(9): 
            lista_txt.append(fila[i].get())
        if lista_txt[0] != '':
            tipoe = 1 if lista_txt[7] == '' else 0 # 0=cuadr_8, 1=triang_6
            nn= 6 if tipoe else 8
            lista_num=[]
            try:
                for i in range(nn+1): 
                    lista_num.append(int(lista_txt[i]))
            except (TypeError, ValueError):
                texto= 'Hay un error en los elems[]. Compruebe'
                texto+=' que cada uno tiene 6 u 8 nodos, que'
                texto+=' las casillas en blanco lo están, etc.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v0)
                return(1)
        e0.tipoe= tipoe
        for i in range(nn):
            j=lista_num[i+1] # es un nodo del elemento
            e0.nodse.append(j)
            e0.xi.append(nodos[j].x)
            e0.yi.append(nodos[j].y)
        elems[lista_num[0]]=e0   
        nelems = len(elems)   # nuevamente actualizar a cada paso
    for ie, e0 in elems.items():
        for i in e0.nodse:
            nodos[i].elemn.append(ie) # informacion redundante en los nodos

    ctes=[]
    try:
        ctes.append(float(entry_young.get()))
        ctes.append(float(entry_poiss.get()))
        ctes.append(float(entry_dilat.get()))
        ctes.append(int(tens_plana.get()))
    except (TypeError, ValueError):
        texto= 'Hay un error en las constantes elásticas.'
        texto+=' Compruebe los valores suministrados.'
        texto+='\n             '+ str(sys.exc_info()[0])
        messagebox.showinfo(message=texto,title='MAL',parent=v0)
        return(1)

    ffx, ffy = [], [] 
    # lista func para ux & lista de func para uy
    for kk in filas_funcs:
        fx,fy= kk[1].get(), kk[2].get()
        if fx != '':
            ffx.append(fx)
            ffy.append(fy)
            nff=len(ffx) 
    if valida_ff(ffx) or valida_ff(ffy): return(1) # validamos por si acaso
            
    cargas=Carga(contorno=[],puntual=[],volumen=[],termica=[])
    
    cargas.contorno.append(0)
    for fila in filas_xraya:
        if fila[1].get() != '':
            nt=int(fila[0].get())  # fila[0] es una StringVar asociada al checkbutton
            try:
                n1,n2,n3= int(fila[1].get()),int(fila[2].get()),int(fila[3].get())
                x1,x2,x3= float(fila[4].get()),float(fila[5].get()),float(fila[6].get())
                y1,y2,y3= float(fila[7].get()),float(fila[8].get()),float(fila[9].get())
            except(TypeError, ValueError):
                texto= 'Hay un error en las cargas de contorno.'
                texto+=' Compruebe los datos suministrados.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v0)
                return(1)
            cargas.contorno.append([nt,n1,n2,n3,x1,x2,x3,y1,y2,y3])
            cargas.contorno[0]=len(cargas.contorno)-1

    cargas.puntual.append(0)
    for fila in filas_fuerzas:
        if fila[0].get() != '':
            try:
                n = int(fila[0].get())
                fuerx, fuery = float(fila[1].get()), float(fila[2].get())
            except(TypeError, ValueError):
                texto= 'Hay un error en las cargas puntuales.'
                texto+=' Compruebe los datos suministrados.'
                texto+='\n             '+ str(sys.exc_info()[0])
                messagebox.showinfo(message=texto,title='MAL',parent=v0)
                return(1)
            cargas.puntual.append([n, fuerx, fuery])
            cargas.puntual[0]=len(cargas.puntual)-1
    
    cargas.volumen=[]
    if filas_xgrande[0].get() != '':
        cargas.volumen.append(filas_xgrande[0].get())
        cargas.volumen.append(filas_xgrande[1].get())
        if valida_ff(cargas.volumen): return(1)

    cargas.termica=[]
    if filas_temperatura[0].get() != '':
        cargas.termica.append( filas_temperatura[0].get())
        if valida_ff(cargas.termica): return(1)
    
    # un extra: asigno variables usadas en subdividir y quiza en dibujar
    extremosx, extremosy = [0.,0.], [0.,0.] # cada uno [min, max]
    for n in nodos.values():
        if n.x < extremosx[0]: extremosx[0]=n.x
        if n.x > extremosx[1]: extremosx[1]=n.x
        if n.y < extremosy[0]: extremosy[0]=n.y
        if n.y > extremosy[1]: extremosy[1]=n.y
    anchox = extremosx[1]-extremosx[0]
    anchoy = extremosy[1]-extremosy[0]


def trocear(elems):
    # divide cada elemento o trozo, en cuatro trozos. No usa los nodos, solo
    # las coord en  e.xi  e.yi
    # Devuelve un diccionario nuevo de trozos (elem_new) en los que e.nodse
    # esta vacio. No modifica los elems de entrada, el programa de llamada
    # lo hara si procede. Los nodos[] globales no se cambian.
    # Notese que elems{} es local aquí
    # Los puntos de trozos contiguos coincidiran -> hay continuidad ok.
    
    ie_current=0
    elems_new= {}  # son los trozos nuevos

    for ie, e in elems.items():
        xi, yi = [],[]
        if e.tipoe == 0:
            xnew, ynew = [], []
            a=[[-0.5,-1], [0.5,-1], [-1,-0.5], [0,-0.5], [1,-0.5], [-0.5,0], [0,0], [0.5,0], [-1,0.5], [0,0.5], [1,0.5], [-0.5,1], [0.5,1] ]
            # son los psi eta de los nodos nuevos (intermedios) a crear
            for p in a:
                x,y=e.dimexy(p[0],p[1])
                xnew.append(x)
                ynew.append(y)

            b= [e.xi[0], xnew[0], e.xi[1], xnew[3], xnew[6], xnew[5], e.xi[7], xnew[2] ]
            c= [e.yi[0], ynew[0], e.yi[1], ynew[3], ynew[6], ynew[5], e.yi[7], ynew[2] ]
            ie_current += 1 # numeramos los elem desde 1
            elems_new[ie_current]=Elem(tipoe=0, xi=b, yi=c)
            
            b= [e.xi[1], xnew[1], e.xi[2], xnew[4], e.xi[3], xnew[7], xnew[6], xnew[3] ]
            c= [e.yi[1], ynew[1], e.yi[2], ynew[4], e.yi[3], ynew[7], ynew[6], ynew[3] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=0, xi=b, yi=c)

            b= [xnew[6], xnew[7], e.xi[3], xnew[10], e.xi[4], xnew[12], e.xi[5], xnew[9] ]
            c= [ynew[6], ynew[7], e.yi[3], ynew[10], e.yi[4], ynew[12], e.yi[5], ynew[9] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=0, xi=b, yi=c)
            
            b= [e.xi[7], xnew[5], xnew[6], xnew[9], e.xi[5], xnew[11], e.xi[6], xnew[8] ]
            c= [e.yi[7], ynew[5], ynew[6], ynew[9], e.yi[5], ynew[11], e.yi[6], ynew[8] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=0, xi=b, yi=c)

        elif e.tipoe ==1:
            xnew, ynew = [], []
            r3= np.sqrt(3.)
            a=[[-0.25,0.25*r3], [0.25,0.25*r3], [0,r3/2], [-0.75, 0.75*r3], [-0.25, 0.75*r3], [0.25, 0.75*r3], [0.75, 0.75*r3], [-0.5, r3], [0.5, r3]]
            # son los psi,eta de los nodos nuevos (intermedios) a crear
            for p in a:
                x,y=e.dimexy(p[0],p[1])
                xnew.append(x)
                ynew.append(y)

            b=[e.xi[0], xnew[1], e.xi[1], xnew[2], e.xi[5], xnew[0] ]
            c=[e.yi[0], ynew[1], e.yi[1], ynew[2], e.yi[5], ynew[0] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=1, xi=b, yi=c)
            
            b=[e.xi[1], xnew[6], e.xi[2], xnew[8], e.xi[3], xnew[5] ]
            c=[e.yi[1], ynew[6], e.yi[2], ynew[8], e.yi[3], ynew[5] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=1, xi=b, yi=c)
            
            b=[e.xi[5], xnew[2], e.xi[1], xnew[5], e.xi[3], xnew[4] ]
            c=[e.yi[5], ynew[2], e.yi[1], ynew[5], e.yi[3], ynew[4] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=1, xi=b, yi=c)
            
            b=[e.xi[5], xnew[4], e.xi[3], xnew[7], e.xi[4], xnew[3] ]
            c=[e.yi[5], ynew[4], e.yi[3], ynew[7], e.yi[4], ynew[3] ]
            ie_current += 1
            elems_new[ie_current]=Elem(tipoe=1, xi=b, yi=c)

    return (elems_new)


def trata_p(carga):
    x,y= nodos[carga[0]].x, nodos[carga[0]].y
    F = np.array([carga[1], carga[2]])  # espero las concentradas en x-y, (NO en n-t) 
    aporta= np.matmul( np.transpose(N_xy(x,y)), F)
    return (aporta)

def trata_s(carga):
    # trato la geometria lineal de 3 ptos como si no fuesen nodos[]
    def H1i(psi):
        coefH1= np.array([[ 1., -3., 2.], [ 0., 4., -4.], [ 0., -1., 2.]])
            # seria mas optimo pero menos claro asignar coefH1 desde fuera
        a= np.matmul( coefH1, np.array([1,psi,psi*psi])  )
        return(a)

    def H1i_psi(psi):   # las derivadas   d H / d psi
        a= np.array([-3.+4.*psi, 4.-8.*psi, -1.+4.*psi])
        return (a)

    def n_ext(psi, x_elem, y_elem):     # vector normal exterior & vector tg
        a= np.matmul( x_elem, H1i_psi(psi) )
        b= np.matmul( y_elem, H1i_psi(psi) )
        c= np.sqrt(a*a+b*b)
        t= np.array([a/c, b/c])
        n= np.array([t[1], -t[0]])
        return(n,t)
    
    a= np.sqrt(5/6.)/6.  # ponemos psi & w de gauss
    w_gaus= np.array([0.5-a, 0.5+a, 0.5+a, 0.5-a]) /2.
    a=2*np.sqrt(6/5.)
    b=np.sqrt( (3+a)/7.)
    a=np.sqrt( (3-a)/7.)
    psi_gaus=[0.5-b/2, 0.5-a/2, 0.5+a/2, 0.5+b/2]

    
    tipoc=carga[0]      # si es nt: 1, si es xy: 0
    xi = np.array([ nodos[carga[1]].x, nodos[carga[2]].x, nodos[carga[3]].x ])
    yi = np.array([ nodos[carga[1]].y, nodos[carga[2]].y, nodos[carga[3]].y ])
    pxi = np.array([carga[4], carga[5], carga[6]])
    pyi = np.array([carga[7], carga[8], carga[9]])
    if tipoc : # si las px py son en realidad pn pt, hay que transformarlas
        a,b, psi = [],[], [0., 0.5, 1.]
        for i in range(3):
            n,t =n_ext(psi[i], xi, yi)
            a.append (pxi[i]*n[0] + pyi[i]*t[0])
            b.append (pxi[i]*n[1] + pyi[i]*t[1])
        pxi, pyi = a, b
    
    aporta=np.zeros(2*nff)    # hago la integral 
    for i in range(4):        # los 4 ptos de Gauss previstos
        x,y = np.matmul( H1i(psi_gaus[i]), xi), np.matmul( H1i(psi_gaus[i]), yi)
        Xraya=np.array([np.matmul( H1i(psi_gaus[i]), pxi), np.matmul(H1i(psi_gaus[i]), pyi)])
        a,b = np.matmul( H1i_psi(psi_gaus[i]), xi), np.matmul( H1i_psi(psi_gaus[i]), yi)
        jaco=np.sqrt(a*a+b*b)
        aporta += np.matmul(np.transpose(N_xy(x,y)), Xraya) * jaco * w_gaus[i]
    
    return (aporta)


def pinta_base():
    # pinta nodos numerados y elementos (del usuario). 
    if  figuras_bien:  plt.axis('equal') # opcion de figuras proporcionadas
    
    for e in elems.values(): e.pinta_elem(colorin='grey', ancholin=1.4, tipolin='-')
    for i,n in nodos.items():
        x,y = n.x, n.y
        plt.plot(x, y, 'kd', markersize=3)
        a = (np.random.random()-0.5) * anchox/25
        b = (np.random.random()-0.5) * anchoy/25
        plt.text(x+a, y+b, str(i), fontsize=9)


def comprobar(elems):
    # salida de texto
    print('#'*47)
    print ('#########   SALIDA DE COMPROBACION   ##########')
    print ('#'*47)
    print('\niNodo     coor x       coor y ')
    for inodo, nodo in nodos.items():
        print('{:4d} {:12.7f} {:12.7f} '.format(inodo, nodo.x, nodo.y))

    for ielem, elem in elems.items():
        print('\niElem: {:3d}   tipo: {:}   '.format(ielem, elem.tipoe))
        print('      inodo     [xy del nodo]:')
        for i in range(n_nodos_tipoe[elem.tipoe]):
            print('     {:4d} {:10.5f} {:10.5f}'.format(
                            elem.nodse[i],elem.xi[i],elem.yi[i]))
    
    print('\nFunciones de aproximación: ', nff)
    lenmedia=0
    for kk in ffx: lenmedia += len(kk)
    lenmedia = int(lenmedia/nff)
    texto= '   i   '+ ' '*(lenmedia//2)+'Para ux'+ ' '*(lenmedia)+ 'Para uy'
    print(texto)
    for i in range(nff):
        print(' {:4d}       {:}       {:}'.format(i, ffx[i], ffy[i]))
    
    
    print ('\n     E_elas       nu_elas      alfa_T      TP?:')
    print (' {:12.4e}  {:12.4e} {:12.4e}    {:}'.format( 
                    ctes[0], ctes[1], ctes[2], ctes[3]))
    
    
    print('\nCargas de contorno (Xraya):')
    if cargas.contorno[0]:
        for i in range(cargas.contorno[0]): print(cargas.contorno[i+1])
    else:
        print('No hay.')

    print('\nCargas puntuales:')
    if cargas.puntual[0]:
        for i in range(cargas.puntual[0]): print(cargas.puntual[i+1])
    else:
        print('No hay.')
    
    print('\nCarga de volumen (Xgrande):')
    if cargas.volumen: 
        print(cargas.volumen)
    else:
        print('No hay.')
            
    print('\nCarga termica (Temperatura escalar):')
    if cargas.termica: 
        print(cargas.termica)
    else:
        print('No hay.')
    
    print ('\n'+'#'*51)
    print ('#########   FIN SALIDA DE COMPROBACION   ##########')
    print ('#'*51 + '\n\n ')


def pinta_comprobar():
    plt.figure('Geometría proporcionada')
    obten_datos()
    # dibujo los nodos:
    for i,n in nodos.items():
        x,y = n.x, n.y
        plt.plot(x, y, 'yd', markersize=3)
        a = (np.random.random()-0.5) * anchox/25
        b = (np.random.random()-0.5) * anchoy/25
        plt.text(x+a, y+b, str(i), fontsize=9)

    # dibujo los bordes en detalle y numero los elem:
    for i,e in elems.items():
        e.pinta_elem(colorin='black', ancholin=0.8, tipolin='')
        x,y= (e.xi[1]+e.xi[5])/2 , (e.yi[1]+e.yi[5])/2
        plt.text(x,y,str(i),fontsize=10, bbox=dict(facecolor='g', 
                edgecolor='white', alpha=0.2, boxstyle='round'))

    plt.show()
    return()


def ver_ff():
    global sin_gui
    quiero_agrupadas= messagebox.askyesno(
        title='¿Agrupar las funciones?',
        message='¿Quiere agrupar las funciones de aproximación en solo dos ventanas?',
        detail='(puede que esas ventanas no quepan en su monitor; si tiene más de 12 funciones por cada componente, quizá prefiera no agruparlas)',
        default='yes', parent=v0 )

    if not sin_gui: obten_datos()
    
    x= np.linspace(extremosx[0], extremosx[1], 99)
    y= np.linspace(extremosy[0], extremosy[1], 99)
    xi, yi = np.meshgrid(x,y)
    x,y=sp.symbols('x,y')
    
    if quiero_agrupadas:
        
        if nff > 3:
            filas_grid=3
            sobrante= nff - 3*(nff//3)
            columnas_grid = nff//3 +1 if sobrante else nff//3
            fig, ejes = plt.subplots(nrows=filas_grid, ncols=columnas_grid, 
                    num='Funciones para ux')

            for i in range(nff):
                ejef, ejec = i-3*(i//3), i//3
                eje=ejes[ejef, ejec]
                
                ff= sp.lambdify([x,y], sp.sympify(ffx[i]), 'numpy')
                zi= ff(xi,yi)
                contornos=eje.contour(xi,yi,zi, levels=[0], colors='red', alpha=0.7)
                eje.clabel(contornos, inline=True, fontsize=8)
                idc=eje.contourf(xi,yi,zi, 34, cmap='Blues', alpha=0.5)
                fig.colorbar(idc, label=ffx[i], ax=eje)
                plt.sca(eje) # hace de este el current axes <= e.pintae() usa plt.plot
                for e in elems.values(): e.pinta_elem(colorin='grey', 
                            ancholin=1.4, tipolin='-')

            fig, ejes = plt.subplots(nrows=filas_grid, ncols=columnas_grid, 
                    num='Funciones para uy')

            for i in range(nff):
                ejef, ejec = i-3*(i//3), i//3
                eje=ejes[ejef, ejec]
                
                ff= sp.lambdify([x,y], sp.sympify(ffy[i]), 'numpy')
                zi= ff(xi,yi)
                contornos=eje.contour(xi,yi,zi, levels=[0], colors='red', alpha=0.7)
                eje.clabel(contornos, inline=True, fontsize=8)
                idc=eje.contourf(xi,yi,zi, 34, cmap='Blues', alpha=0.5)
                fig.colorbar(idc, label=ffy[i], ax=eje)
                plt.sca(eje) # hace de este el current axes <= e.pintae() usa plt.plot
                for e in elems.values(): e.pinta_elem(colorin='grey', 
                            ancholin=1.4, tipolin='-')
        else:
            fig, ejes = plt.subplots(nrows=nff, ncols=1, num='Funciones para ux')
            
            for i in range(nff):
                eje=ejes[i]
                ff= sp.lambdify([x,y], sp.sympify(ffx[i]), 'numpy')
                zi= ff(xi,yi)
                contornos=eje.contour(xi,yi,zi, levels=[0], colors='red', alpha=0.7)
                eje.clabel(contornos, inline=True, fontsize=8)
                idc=eje.contourf(xi,yi,zi, 34, cmap='Blues', alpha=0.5)
                fig.colorbar(idc, label=ffx[i], ax=eje)
                plt.sca(eje) # hace de este el current axes <= e.pintae() usa plt.plot
                for e in elems.values(): e.pinta_elem(colorin='grey', 
                            ancholin=1.4, tipolin='-')
                            
            fig, ejes = plt.subplots(nrows=nff, ncols=1, num='Funciones para uy')
            
            for i in range(nff):
                eje=ejes[i]
                ff= sp.lambdify([x,y], sp.sympify(ffy[i]), 'numpy')
                zi= ff(xi,yi)
                contornos=eje.contour(xi,yi,zi, levels=[0], colors='red', alpha=0.7)
                eje.clabel(contornos, inline=True, fontsize=8)
                idc=eje.contourf(xi,yi,zi, 34, cmap='Blues', alpha=0.5)
                fig.colorbar(idc, label=ffy[i], ax=eje)
                plt.sca(eje) # hace de este el current axes <= e.pintae() usa plt.plot
                for e in elems.values(): e.pinta_elem(colorin='grey', 
                            ancholin=1.4, tipolin='-')

    else:  # las quiero sin agrupar

        for i in range(nff):
            name= 'ux'+str(i)
            plt.figure(name)
            ff= sp.lambdify([x,y], sp.sympify(ffx[i]), 'numpy')
            zi= ff(xi,yi)
            contornos=plt.contour(xi,yi,zi, levels=[0], colors='red', alpha=0.7)
            plt.clabel(contornos, inline=True, fontsize=8)
            plt.contourf(xi,yi,zi, 34, cmap='Blues', alpha=0.5)
            plt.colorbar(label=ffx[i])
            pinta_base()
            
            name= 'uy'+str(i)
            plt.figure(name)
            ff= sp.lambdify([x,y], sp.sympify(ffy[i]), 'numpy')
            zi= ff(xi,yi)
            contornos=plt.contour(xi,yi,zi, levels=[0], colors='red', alpha=0.7)
            plt.clabel(contornos, inline=True, fontsize=8)
            plt.contourf(xi,yi,zi, 34, cmap='Blues', alpha=0.5)
            plt.colorbar(label=ffy[i])
            pinta_base()

    if not sin_gui: plt.show() # si estamos sin_gui que siga de corrido



def a_guardar():
    obten_datos()
    #nfcompleto = filedialog.asksaveasfilename(parent=v0,initialdir=dir_home,
    #    title='Guardar problema')
    nfcompleto = filedialog.asksaveasfilename(parent=v0, title='Guardar problema')
    if bool(nfcompleto):
        n_f = path.basename(nfcompleto)
        v0.title('Galerkin_Appr v0.7 - '+n_f)

        f=open(nfcompleto,'w')
        
        print('Identificación: ', nfcompleto,file=f)
        print('\nNum_de_nodos     Num_de_elementos', file=f)
        print(' {:8d}  {:15d}'.format(nnodos, nelems),file=f)
        
        print('\niNodo      coor_x      coor_y', file=f)
        for i,n in nodos.items():
            print('{:3d}  {:10.5f}  {:10.5f}'.format(i, n.x, n.y),file=f)
        print('\niElem     tipoE            nodos', file=f)
        for i,e in elems.items():
            s = '  ' + str(i) + '   ' + str(e.tipoe) + '     '
            for n in e.nodse:
                s += '  ' + str(n)
            print(s, file=f)
        print('\nFunciones de aproximacion:', file=f)
        print ('  {:3d} '.format(len(ffx)), file=f)
        print('  Para ux:', file=f)
        for s in ffx:
            print('    '+s, file=f)
        print('  Para uy:', file=f)
        for s in ffy:
            print('    '+s, file=f)
        print('\n     E        nu      alfa_temp     TP/DP', file=f)
        a='TP' if ctes[3] else 'DP'
        print('{:8.3f}  {:8.3f}  {:8.3f}          {:8}'.format(
            ctes[0], ctes[1], ctes[2], a), file=f)
    
        print('\nCargas: contorno -> puntual -> volumen(1ó0) -> termica(1ó0)', file=f)

        print('  contorno {:3d}'.format(len(cargas.contorno)-1), file=f)        
        if cargas.contorno:
            for i in range(len(cargas.contorno)-1) :
                s='    nt' if cargas.contorno[i+1][0] else '    xy'
                for j in cargas.contorno[i+1][1:]:
                    s += ' ' + str(j)
                print(s, file=f)

        print('  puntual  {:3d}'.format(len(cargas.puntual)-1), file=f)
        if cargas.puntual:
            for i in range(len(cargas.puntual)-1) :
                print('  {:4d}  {:8.3f}  {:8.3f}'.format(cargas.puntual[i+1][0],
                cargas.puntual[i+1][1],cargas.puntual[i+1][2]), file=f)
        
        print('  volumen  {:3d}'.format(int(len(cargas.volumen)/2)), file=f)
        if cargas.volumen:
            print('  ' + cargas.volumen[0], file=f)
            print('  ' + cargas.volumen[1], file=f)

        print('  termica  {:3d}'.format(int(len(cargas.termica))), file=f)
        if cargas.termica:
            print('  ' + cargas.termica[0], file=f)
        
        print('# ___fin___', file=f)
        f.close()
        
        
    
def a_cargar():
    ####  Lectura de fichero  ####
    global elige, n_f, nfcompleto, nodos, elems, nnodos, nelems
    global ctes, nff, ffx, ffy, cargas, sin_gui, niveles_sigma_bak
    nodos, elems, ctes, ffx, ffy, cargas = {}, {}, [], [], [], Carga([],[],[],[])
    for i in range(7):
        niveles_sigma_bak[i][0], niveles_sigma_bak[i][1] = '', ''

    #nfcompleto=filedialog.askopenfilename(parent=v0, initialdir=dir_home,
    #            title='Abrir archivo')
    nfcompleto=filedialog.askopenfilename(parent=v0, title='Abrir archivo')
    if not nfcompleto: return()
    n_f= path.basename(nfcompleto)
    v0.title('Galerkin Appr - '+n_f)

    f=open(nfcompleto, 'r')
    
    ### OJO OJO: guardo en las variables fundamentales (nodos, elems...) ###
    ### los contenidos del fichero para no liarme, pero esto no          ###
    ### sobrevivira porque calcula() llamara a obten_datos(), Asi que lo ###
    ### importante es escribir esto bien en el GUI, como hago luego.     ###

    kk= f.readline()    # Titulo o lo que quieras
    kk= f.readline()    #linea en blanco
    kk= f.readline()    # texto Num nodos, Num elems
    kk= f.readline().split()
    nnodos, nelems = int(kk[0]), int(kk[1])
    kk= f.readline()    #linea en blanco

    kk= f.readline()    # texto iNodo, coor x, coor y
    for i in range(nnodos):
        kk=f.readline().split()
        j=int(kk[0])
        nodos[j]= Nodo(0.,0.,0.,0.,[],[])
        nodos[j].x, nodos[j].y =   float(kk[1]), float(kk[2])
        nodos[j].elemn =[]  # aprovechamos para inicializar n.elemn
    kk= f.readline()    #linea en blanco


    kk=f.readline()       # texto iElem, tipoE, nodos del elemento
    for i in range(nelems):
        kk=f.readline().split()
        j=int(kk[0])
        elems[j]= Elem()
        elems[j].tipoe = int(kk[1])
        elems[j].nodse = kk[2:]
        for k in range(n_nodos_tipoe[elems[j].tipoe]):
            elems[j].nodse[k]= int(elems[j].nodse[k])

    for ielem,e in elems.items():  # rellenar n.elemn[]
        for inode in e.nodse:
            nodos[inode].elemn.append(ielem) # a ver si funciona sin inicializar aqui
        if sin_gui: # si estamos sin_gui, rellenar e.xi, e.yi aqui (no hay obten_datos):
            e.xi, e.yi = [], []
            for inode in e.nodse:
                e.xi.append(nodos[inode].x)
                e.yi.append(nodos[inode].y)

    kk= f.readline() 
    kk= f.readline()    # texto Funciones de aproximacion nff - ffx - ffy
    nff= int( f.readline()) 
    ffx, ffy = [], []

    kk=f.readline()     # texto 'para ux'
    for i in range(nff):
        ffx.append(f.readline().strip())  # para que no lea el CR ni espacios!
    kk= f.readline()    # texto 'para uy'
    for i in range(nff):
        ffy.append(f.readline().strip())

    kk= f.readline() 
    kk= f.readline()    # texto 'constantes elasticas E, nu, alfa_dilat, TP/DP
    kk= f.readline().split()
    a = 1 if kk[3]=='TP' else 0
    ctes=[float(kk[0]),float(kk[1]), float(kk[2]), a]

    kk= f.readline() 
    kk= f.readline()    # texto 'cargas cargas por orden: distr, punt, volum, termicas'

    cargas.contorno = []
    n_cargas_s = int( f.readline().split()[1])  # numero de cargas s (de contorno)
    cargas.contorno.append(n_cargas_s)
    for i in range(n_cargas_s):
        a = f.readline().split()
        local =1 if a[0]=='nt' else 0
        b= [ local, int(a[1]), int(a[2]), int(a[3]), float(a[4]), float(a[5]), 
                           float(a[6]), float(a[7]), float(a[8]), float(a[9])]
        cargas.contorno.append(b)

    cargas.puntual = []
    n_cargas_p = int( f.readline().split()[1])  # numero de cargas p (puntuales)
    cargas.puntual.append(n_cargas_p)
    for i in range(n_cargas_p):
        a = f.readline().split()  # inodo, fuerza_x, fuerza_y
        b= [ int(a[0]), float(a[1]), float(a[2]) ]
        cargas.puntual.append(b)

    # habra una o ninguna carga de V, en dos componentes 
    cargas.volumen=[]
    n_cargas_V = int( f.readline().split()[1])  # numero de cargas de vol (0 o 1)
    if n_cargas_V:
        cargas.volumen.append(f.readline().strip() )
        cargas.volumen.append(f.readline().strip() ) # lee X_x, X_y, como funciones (x,y)

    # habra una o ninguna carga termica; lee T(xy), hay que componer [eps_0] aparte
    cargas.termica=[]
    n_cargas_T = int( f.readline().split()[1])  # numero de cargas T (0 o 1)
    if n_cargas_T: cargas.termica.append(f.readline().strip() )

    f.close()
    
    if sin_gui: return()  # es que estamos obviando gui 
    
    #### Fin de lecturas.        ####
    #### Ahora rellenamos el GUI ####


    a_borrar()
    iposi=0
    for i,n  in nodos.items():
        filas_nodos[iposi][0].insert(0,str(i))
        filas_nodos[iposi][1].insert(0,str(n.x))
        filas_nodos[iposi][2].insert(0,str(n.y))
        iposi +=1

    # ielem, nodos  (el tipoe es una cosa interna)
    iposi=0    
    for i,e in elems.items():
        filas_elems[iposi][0].delete(0,'end')
        filas_elems[iposi][0].insert(0,str(i))
        for j in range(n_nodos_tipoe[e.tipoe]):
            filas_elems[iposi][j+1].delete(0,'end')
            filas_elems[iposi][j+1].insert(0,str(e.nodse[j]))
        iposi +=1
    
    # ctes elasticas
    entry_young.delete(0,'end')
    entry_poiss.delete(0,'end')
    entry_dilat.delete(0,'end')
    entry_young.insert(0,str(ctes[0]))
    entry_poiss.insert(0,str(ctes[1]))
    entry_dilat.insert(0,str(ctes[2]))
    tens_plana.set(str(ctes[3]))

    # funciones de aproximacion 
    for i in range(nff):
        filas_funcs[i][1].delete(0,'end')
        filas_funcs[i][2].delete(0,'end')
        filas_funcs[i][1].insert(0, ffx[i]) # para ux
        filas_funcs[i][2].insert(0, ffy[i]) # para uy



    # cargas
    for i in range(cargas.contorno[0]):
        filas_xraya[i][0].set(str(cargas.contorno[i+1][0]))
        
        filas_xraya[i][1].insert(0,str(cargas.contorno[i+1][1]))
        filas_xraya[i][2].insert(0,str(cargas.contorno[i+1][2]))
        filas_xraya[i][3].insert(0,str(cargas.contorno[i+1][3]))
        filas_xraya[i][4].insert(0,str(cargas.contorno[i+1][4]))
        filas_xraya[i][5].insert(0,str(cargas.contorno[i+1][5]))
        filas_xraya[i][6].insert(0,str(cargas.contorno[i+1][6]))
        filas_xraya[i][7].insert(0,str(cargas.contorno[i+1][7]))
        filas_xraya[i][8].insert(0,str(cargas.contorno[i+1][8]))
        filas_xraya[i][9].insert(0,str(cargas.contorno[i+1][9]))

    for i in range(cargas.puntual[0]):
        filas_fuerzas[i][0].insert(0,str(cargas.puntual[i+1][0]))
        filas_fuerzas[i][1].insert(0,str(cargas.puntual[i+1][1]))
        filas_fuerzas[i][2].insert(0,str(cargas.puntual[i+1][2]))

    if cargas.volumen:
        filas_xgrande[0].insert(0,str(cargas.volumen[0]))
        filas_xgrande[1].insert(0,str(cargas.volumen[1]))
    
    if cargas.termica:
        filas_temperatura[0].insert(0,str(cargas.termica[0]))
    
    return()


def borrar_una_lista(lista_de_casillas): # sirve solo si es array 2d de entrys
    for fila in lista_de_casillas:
        for j in range(len(fila)):
            fila[j].delete(0,'end')


def a_limpiar(): # ventana de borrar por partes
    
    def hecholimpiar():
        
        if lista_bool[0].get(): 
            borrar_una_lista(filas_nodos)
            borrar_una_lista(filas_elems)
            entry_young.delete(0,'end')
            entry_poiss.delete(0,'end')
            entry_dilat.delete(0,'end')
            tens_plana.set('1')

        if lista_bool[1].get():
            borrar_una_lista(filas_funcs)

        if lista_bool[2].get():
            borrar_una_lista(filas_fuerzas)
            for fila in filas_xraya: # es especial por tener un pulsador
                fila[0].set('0')
                for i in range(1,10): fila[i].delete(0,'end')
            filas_xgrande[0].delete(0,'end')
            filas_xgrande[1].delete(0,'end')
            filas_temperatura[0].delete(0,'end')
            
        v3.destroy()


    v3=Toplevel(v0)
    v3.title('Limpiar campos')
    v3.geometry('+300+200')

    lista_bool=[BooleanVar(), BooleanVar(), BooleanVar()]
            # son variables de los botones
    
    a= 'Marque los campos de entradas\n'
    a+='   que quiere dejar en blanco:'
    ttk.Label(v3, text=a, background='#EDECEB').grid(
            row=0, column=0, columnspan=2)
            
    ttk.Checkbutton(v3, variable=lista_bool[0]).grid(
            row=1, column=0)
    ttk.Label(v3, text='Dominio', background='#EDECEB').grid(
            row=1, column=1, sticky='w')

    ttk.Checkbutton(v3, variable=lista_bool[1]).grid(
            row=2, column=0)
    ttk.Label(v3, text='Funciones', background='#EDECEB').grid(
            row=2, column=1, sticky='w')
    
    ttk.Checkbutton(v3, variable=lista_bool[2]).grid(
            row=3, column=0)
    ttk.Label(v3, text='Cargas', background='#EDECEB').grid(
            row=3, column=1, sticky='w')
    
    ttk.Button(v3, text='Proceder', command=hecholimpiar).grid(
            row=4, column=0, columnspan=2)

    for hijo in v3.winfo_children(): hijo.grid_configure(padx=6,
        pady=6)

    #v3.geometry('+410-70')
    v3.focus()
    v3.mainloop()


def a_borrar():
    borrar_una_lista(filas_nodos)
    borrar_una_lista(filas_elems)
    entry_young.delete(0,'end')
    entry_poiss.delete(0,'end')
    entry_dilat.delete(0,'end')
    tens_plana.set('1')
    
    borrar_una_lista(filas_funcs)
    
    borrar_una_lista(filas_fuerzas)
    for fila in filas_xraya: # es especial por tener un pulsador
        fila[0].set('0')
        for i in range(1,10): fila[i].delete(0,'end')
    filas_xgrande[0].delete(0,'end')
    filas_xgrande[1].delete(0,'end')
    filas_temperatura[0].delete(0,'end')
    

def a_salir():
    if messagebox.askokcancel(message='¿Quiere cerrar el programa? ',
                detail='Los datos no guardados se perderán.', default='cancel',
                icon='question', title='Confirmacion:',parent=v0) :
        for i in plt.get_fignums(): plt.close(i)
        v0.destroy()
        try: # si es linux que cierre el terminal tambien
            kill(getppid(), signal.SIGHUP)        
        except:
            exit()
        # con extension .pyw no saca terminal.

def cierraplots():
    global v5
    for i in plt.get_fignums(): plt.close(i)
    try:
        v5.destroy() # tambien la ventana de resultados
    except:
        pass



# Ventana de opciones 

def entraOpc():
    global mostrar_integracion, precision_integrar, figuras_bien, sigma_polar
    
    def hechoOpc():
        global mostrar_integracion, precision_integrar, figuras_bien, sigma_polar
        
        mostrar_integracion= int(mostrar_di.get())
        precision_integrar= int(precis_di.get())
        figuras_bien= int(figsb_di.get())
        sigma_polar = int(sigma_polar_di.get())
        
        npgaus= kk1*4**precision_integrar
        if  kk1 > 2350: # deja hacer el ej con precision maxima
            texto = 'Su configuración implica un número\n'
            texto = 'excesivo de puntos de muestro para\n'
            texto = 'la integración ({:4d}). Se trata de\n'.format (npgaus)
            texto+= 'una limitación impuesta por seguridad.\n'
            texto+= 'Por favor elija una integración menos\n'
            texto+= 'exigente.'
            messagebox.showinfo(message=texto)
            return()

        v4.destroy()

    v4=Toplevel(v0)
    v4.title('Opciones')

    mostrar_di= StringVar()
    precis_di = StringVar()
    figsb_di  = StringVar()
    sigma_polar_di  = StringVar()
    
    
    try:
        mostrar_di= StringVar(value=str(mostrar_integracion))
    except (ValueError,NameError):
        mostrar_di= StringVar()
        
    try:
        precis_di = StringVar(value=str(precision_integrar))
    except (ValueError,NameError):
        precis_di = StringVar()
    
    try:
        figsb_di  = StringVar(value=str(figuras_bien))
    except (ValueError,NameError):
        figsb_di  = StringVar()

    try:
        sigma_polar_di  = StringVar(value=str(sigma_polar))
    except (ValueError,NameError):
        sigma_polar_di  = StringVar()
    
    a= '\nMarque esta casilla para obtener una figura \n'
    a+='informativa sobre cómo se ha realizado la\n' 
    a+='integración cada vez que recalcule.\n'
    ttk.Label(v4, text='Mostrar info de integración.',
        background='#EDECEB').grid(row=0, column=1, sticky='w')
    mostrar_bt=ttk.Checkbutton(v4, variable=mostrar_di)
    mostrar_bt.grid(row=0, column=0, sticky='ns')
    CreateToolTip(mostrar_bt,a)
    
    ttk.Separator(v4, orient=HORIZONTAL).grid(
            row=1, column=1,sticky='ew')

    kk1=0
    for e in elems.values():
        a=6 if e.tipoe else 9
        kk1 += a
    kk1 *= 4

    a =  'Precisión de las integraciones: La opción por\n'
    a += 'defecto usará unos {:4d} puntos de muestreo.'.format(kk1)
    ttk.Label(v4, text=a, foreground='green',  background='#EDECEB').grid(
            row=2, column=1, sticky='w')
    ttk.Separator(v4, orient=VERTICAL).grid(
            row=2, column=0,sticky='ns')
    
    a =  '\nLa integración se realiza por cuadratura parabólica.\n'
    a += 'El programa calcula un número de puntos de muestreo\n'
    a += 'adecuado a este problema asumiendo que las funciones\n'
    a += 'de aproximación evolucionan suavemente en el dominio.\n'
    a += 'Si estima que esa opción puede ser insuficiente \n'
    a += 'puede multiplicar ese número por 4, 16 o 64:\n'

    precis0_bt=ttk.Radiobutton(v4, variable=precis_di, value='0')
    precis0_bt.grid(row=3, column=0)
    ttk.Label(v4, text='Precisión por defecto',  background='#EDECEB').grid(
            row=3, column=1, sticky='w')
    CreateToolTip(precis0_bt,a)
            
    precis1_bt=ttk.Radiobutton(v4, variable=precis_di, value='1')
    precis1_bt.grid(row=4, column=0)
    ttk.Label(v4, text='Precisión mejorada (x4)', background='#EDECEB').grid(
            row=4, column=1, sticky='w')
    CreateToolTip(precis1_bt,a)

    precis2_bt=ttk.Radiobutton(v4, variable=precis_di, value='2')
    precis2_bt.grid(row=5, column=0)
    ttk.Label(v4,text='Precisión loca (x16)', background='#EDECEB').grid(
            row=5, column=1, sticky='w')
    CreateToolTip(precis2_bt,a)

    precis3_bt=ttk.Radiobutton(v4, variable=precis_di, value='3')
    precis3_bt.grid(row=6, column=0)
    ttk.Label(v4,text='Precisión muy loca (x64)', background='#EDECEB').grid(
            row=6, column=1, sticky='w')
    CreateToolTip(precis3_bt,a)
    
    ttk.Separator(v4, orient=HORIZONTAL).grid(
            row=7, column=1, sticky='ew')

    a= 'Modo de visualización de las figuras de\n'
    a+='resultados (tensiones y desplazamientos).'
    ttk.Label(v4, text=a, foreground='green',  background='#EDECEB').grid(
            row=8, column=1, sticky='w')
    ttk.Separator(v4, orient=VERTICAL).grid(
            row=8, column=0, sticky='ns')
    ttk.Label(v4,text='Figuras en proporciones reales', background='#EDECEB').grid(
            row=9, column=1, sticky='w')
    figsb_bt=ttk.Checkbutton(v4, variable=figsb_di)
    figsb_bt.grid(row=9, column=0, sticky='ns')
    a= 'Por defecto las figuras de resultados se trazan \n'
    a+='sin respetar sus proporciones ancho/alto para \n'
    a+='aprovechar mejor el espacio de visualización. \n'
    a+='Usted puede si lo desea cambiar manualmente la \n'
    a+='geometría de la ventana, y el dibujo que contiene \n'
    a+='también lo hará. No obstante si prefiere ver la \n'
    a+='geometría con sus proporciones reales desde el \n'
    a+='principio, use esta opción.'
    CreateToolTip(figsb_bt,a)
    
    
    ttk.Separator(v4, orient=HORIZONTAL).grid(
        row=10, column=1,sticky='ew')


    a= 'Trazados de tensiones en coordenadas\n'
    a+='polares (sigma_rr, sigma_tt, sigma_rt)'
    ttk.Label(v4, text=a, foreground='green',  background='#EDECEB').grid(
            row=11, column=1, sticky='w')
    ttk.Separator(v4, orient=VERTICAL).grid(
            row=11, column=0, sticky='ns')
    ttk.Label(v4,text='Obtener tensiones en polares', background='#EDECEB').grid(
            row=12, column=1, sticky='w')
    sigma_polar_bt=ttk.Checkbutton(v4, variable=sigma_polar_di)
    sigma_polar_bt.grid(row=12, column=0, sticky='ns')
    a= 'Si en este problema le interesan las tensiones \n'
    a+='en coordenadas polares, marque esta opción para \n'
    a+='que se calculen. La salida presentará gráficos \n'
    a+='adicionales de estas tensiones.'
    CreateToolTip(sigma_polar_bt, a)
    
    ttk.Separator(v4, orient=HORIZONTAL).grid(
        row=13, column=1,sticky='ew')
    
    ttk.Button(v4, text='Hecho', command=hechoOpc).grid(
            row=14, column=0, columnspan=2)

    for hijo in v4.winfo_children(): hijo.grid_configure(padx=6,
        pady=6)

    v4.geometry('+410-70')
    v4.focus()
    v4.mainloop()


def sigma_vm(s):
    # en s[] las tres componentes xy de tension en el plano
    s33=0 if ctes[3] else ctes[1]*(s[0]+s[1])
    s_vm = (s[0]-s[1])**2+ (s[0]-s33)**2+ (s[1]-s33)**2 + 6*s[2]**2
    s_vm = np.sqrt(s_vm/2)
    return (s_vm)

def a_polares(x,y,sxx,syy,sxy):
    # en  s[] las tres componentes de tension en el punto x,y
    r= np.hypot(x,y)
    c, s = x/r , y/r
    cc, ss, sc = c*c, s*s, s*c
    srr=  cc*sxx + ss*syy + 2*sc*sxy
    stt=  ss*sxx + cc*syy - 2*sc*sxy
    srt= -sc*(sxx-syy) - (ss-cc)*sxy
    return(np.asarray([srr,stt,srt]))

def doscifras(a):
    return(float('{:.2g}'.format(a)))



def salida_grafica():

    global niveles_sigma, niveles_sigma_bak
    
    sigma = lambda x,y: np.matmul( 
                    matrizD, (np.matmul(LN_xy(x,y), a_despl) - eps0_xy(x,y)))

    xi, yi, triangulos, max_pto = np.array([]), np.array([]), np.array([]), 0
    for e in elems.values():
        xi_aport, yi_aport = e.dimexy(e.triangulos_[e.tipoe].x,e.triangulos_[e.tipoe].y)
        xi= np.append(xi,xi_aport, axis=0)
        yi= np.append(yi,yi_aport, axis=0)
        t= e.triangulos_[e.tipoe].triangles + max_pto
        if triangulos.size:
            triangulos= np.append(triangulos, t, axis=0)
        else:
            triangulos = t # inicializa los triangulos la 1a vez
        max_pto += len (xi_aport)

    xi, yi = np.array(xi), np.array(yi)
    triang = matplotlib.tri.Triangulation(xi, yi, triangles=triangulos) 

    #sigmai = sigma(xi,yi)
    sigmai=[[],[],[]] # codigo para arreglar que lo anterior ya no pita- may2024
    for i in range(len(xi)):
        sigma_pto= sigma(xi[i], yi[i])
        sigmai[0].append(sigma_pto[0])
        sigmai[1].append(sigma_pto[1])
        sigmai[2].append(sigma_pto[2])
    sigmai=np.array(sigmai)

    plt.figure('Sigma_xx')
    fig, ax = plt.gcf(), plt.gca()
    if niveles_sigma[0][0] == '' or niveles_sigma[0][1] == '':
        niveles_sigma_bak[0][0] = doscifras(min(sigmai[0]))
        niveles_sigma_bak[0][1] = doscifras(max(sigmai[0]))
        niveles, niveles2 = 21, 5
    else:
        niveles=np.linspace(niveles_sigma[0][0],niveles_sigma[0][1],21) 
        niveles2 = niveles[0::5]
    rellenos= ax.tricontourf(triang, sigmai[0], niveles, cmap='terrain',alpha=0.5)
    cb=fig.colorbar(rellenos, label='sigma_xx')
    contornos=ax.tricontour(triang, sigmai[0], niveles2, colors='red', alpha=0.7)
    contornos.clabel(inline=True, fontsize=8)
    pinta_base()


    plt.figure('Sigma_yy')
    fig, ax = plt.gcf(), plt.gca()
    if niveles_sigma[1][0] == '' or niveles_sigma[1][1] == '':
        niveles_sigma_bak[1][0] = doscifras(min(sigmai[1]))
        niveles_sigma_bak[1][1] = doscifras(max(sigmai[1]))
        niveles, niveles2 = 21, 5
    else:
        niveles=np.linspace(niveles_sigma[1][0],niveles_sigma[1][1],21) 
        niveles2 = niveles[0::5]
    rellenos= ax.tricontourf(triang, sigmai[1], niveles, cmap='terrain',alpha=0.5)
    fig.colorbar(rellenos, label='sigma_yy')
    contornos=ax.tricontour(triang, sigmai[1], niveles2, colors='red', alpha=0.7)
    contornos.clabel(inline=True, fontsize=8)
    pinta_base()

    plt.figure('Sigma_xy')
    fig, ax = plt.gcf(), plt.gca()
    if niveles_sigma[2][0] == '' or niveles_sigma[2][1] == '':
        niveles_sigma_bak[2][0] = doscifras(min(sigmai[2]))
        niveles_sigma_bak[2][1] = doscifras(max(sigmai[2]))
        niveles, niveles2 = 21, 5
    else:
        niveles=np.linspace(niveles_sigma[2][0],niveles_sigma[2][1],21) 
        niveles2 = niveles[0::5]
    rellenos= ax.tricontourf(triang, sigmai[2], niveles, cmap='terrain',alpha=0.5)
    fig.colorbar(rellenos, label='sigma_xy')
    contornos=ax.tricontour(triang, sigmai[2], niveles2, colors='red', alpha=0.7)
    contornos.clabel(inline=True, fontsize=8)
    pinta_base()
    
    plt.figure('Tension de von Mises')
    fig, ax = plt.gcf(), plt.gca()
    sigmai_vm=sigma_vm(sigmai)
    
    if niveles_sigma[3][0] == '' or niveles_sigma[3][1] == '':
        niveles_sigma_bak[3][0] = doscifras(min(sigmai_vm))
        niveles_sigma_bak[3][1] = doscifras(max(sigmai_vm))
        niveles, niveles2 = 21, 5
    else:
        niveles=np.linspace(niveles_sigma[3][0],niveles_sigma[3][1],21)
        niveles2 = niveles[0::5]
    rellenos= ax.tricontourf(triang, sigmai_vm, niveles, cmap='terrain',alpha=0.5)
    fig.colorbar(rellenos, label='sigma_vM')
    contornos=ax.tricontour(triang, sigmai_vm, niveles2, colors='red', alpha=0.7)
    contornos.clabel(inline=True, fontsize=8)
    pinta_base()
    
    if sigma_polar:
        sigmai_pol= a_polares(xi, yi, sigmai[0], sigmai[1], sigmai[2])

        plt.figure('Sigma_rr')
        fig, ax = plt.gcf(), plt.gca()
        if niveles_sigma[4][0] == '' or niveles_sigma[4][1] == '':
            niveles_sigma_bak[4][0] = doscifras(min(sigmai_pol[0]))
            niveles_sigma_bak[4][1] = doscifras(max(sigmai_pol[0]))
            niveles, niveles2 = 21, 5
        else:
            niveles=np.linspace(niveles_sigma[4][0],niveles_sigma[4][1],21) 
            niveles2 = niveles[0::5]
        rellenos= ax.tricontourf(triang, sigmai_pol[0], niveles, cmap='terrain',alpha=0.5)
        fig.colorbar(rellenos, label='sigma_rr')
        contornos=ax.tricontour(triang, sigmai_pol[0], niveles2, colors='red', alpha=0.7)
        contornos.clabel(inline=True, fontsize=8)
        pinta_base()

        plt.figure('Sigma_tt')
        fig, ax = plt.gcf(), plt.gca()

        if niveles_sigma[5][0] == '' or niveles_sigma[5][1] == '':
            niveles_sigma_bak[5][0] = doscifras(min(sigmai_pol[1]))
            niveles_sigma_bak[5][1] = doscifras(max(sigmai_pol[1]))
            niveles, niveles2 = 21, 5
        else:
            niveles=np.linspace(niveles_sigma[5][0],niveles_sigma[5][1],21) 
            niveles2 = niveles[0::5]
        rellenos= ax.tricontourf(triang, sigmai_pol[1], niveles, cmap='terrain',alpha=0.5)
        fig.colorbar(rellenos, label='sigma_tt')
        contornos=ax.tricontour(triang, sigmai_pol[1], niveles2, colors='red', alpha=0.7)
        contornos.clabel(inline=True, fontsize=8)
        pinta_base()


        plt.figure('Sigma_rt')
        fig, ax = plt.gcf(), plt.gca()

        if niveles_sigma[6][0] == '' or niveles_sigma[6][1] == '':
            niveles_sigma_bak[6][0] = doscifras(min(sigmai_pol[2]))
            niveles_sigma_bak[6][1] = doscifras(max(sigmai_pol[2]))
            niveles, niveles2 = 21, 5
        else:
            niveles=np.linspace(niveles_sigma[6][0],niveles_sigma[6][1],21) 
            niveles2 = niveles[0::5]
        rellenos= ax.tricontourf(triang, sigmai_pol[2], niveles, cmap='terrain',alpha=0.5)
        fig.colorbar(rellenos, label='sigma_rt')
        contornos=ax.tricontour(triang, sigmai_pol[2], niveles2, colors='red', alpha=0.7)
        contornos.clabel(inline=True, fontsize=8)
        pinta_base()
    
    #plt.figure('Triángulos para dibujar')
    #plt.triplot(triang, 'k.-', lw=0.2)
    #pinta_base()
    

    # campo de desplazamientos
    
    def esta_dentro(x,y):
        for i,e in elems.items():
            v2x, v2y = e.xi-x, e.yi-y
            v1x, v1y = np.roll(v2x,1), np.roll(v2y,1)
            p=np.cross([v1x, v1y], [v2x,v2y], axisa=0, axisb=0 )
            area= np.sum(p)
            if all (a_sect>0 for a_sect in p):
                return (True)
            else: # dejamos algo de margen, si no en lados curvos quita demasiado
                area_mal= 0.
                for a in p:
                    if a<0. : area_mal += a
                if abs(area_mal/area) < 0.01:
                    return(True)
        return(False)
        
    mx,my= anchox/18., anchoy/18. # mas margen, sitio para la deformada
    x= np.linspace(extremosx[0]-mx, extremosx[1]+mx, 16)
    y= np.linspace(extremosy[0]-my, extremosy[1]+my, 16)
    xi, yi = [] , []
    for kk1 in x:   # meshgrid me lia ahora, quiver no admite mask
        for kk2 in y:
            if esta_dentro (kk1,kk2):
                xi.append(kk1)
                yi.append(kk2)
    
    xi= np.asarray(xi)
    yi= np.asarray(yi)

    u_despl= lambda x,y: np.matmul(N_xy(x,y), a_despl)
    #u_despli= u_despl(xi,yi)
    u_despli=[[],[]] # como antes, para arreglar que lo anterior ya no va may2024
    for i in range(len(xi)):
        u_pto= u_despl(xi[i], yi[i])
        u_despli[0].append(u_pto[0])
        u_despli[1].append(u_pto[1])
    u_despli= np.array(u_despli)

    plt.figure('Desplazamientos')
    plt.quiver(xi,yi, u_despli[0], u_despli[1], units='dots', width=1, color='b') 
            # eso son las flechas
    
    umax = max ( abs(np.amax(u_despli[0])), abs(np.amin(u_despli[0])) ) # para escala deformada
    uscal= 0.02*(anchox+anchoy)/umax
    for e in elems.values(): 
        e.pinta_elem(colorin='m', tipolin='-', ancholin=1, desplacin=True, uscal=uscal)
    pinta_base()

    
    plt.show(block=False)
    return()


def calcula():
    global nnodos,nelems,nff, nodos,elems,ffx,ffy, ctes, cargas
    global matrizD, N_xy, LN_xy, a_despl, eps0_xy 
    global extremosx, extremosy, anchox, anchoy, v5
    
    if not sin_gui: # si estamos con gui que obtenga datos de ahi
        hay_err = obten_datos()
        if hay_err: return(1)
    
    def haz_D():
        # devuelve la matriz D. No hace nu_fake, usa E & nu originales
        E, nu= ctes[0], ctes[1]
        D= np.array
        if ctes[3]: # sera =1, TP
            a=E/(1-nu*nu)
            D = np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
            D *= a
        else: # sera =0, DP
            a, b = E/((1.+nu)*(1.-2*nu)), (1.-nu)
            D = np.array([[b,nu,0],[nu,b,0],[0,0,(1-2*nu)/2]])
            D *= a
        return(D)

    def integrar_en_V():
        # subdivision e integracion por cuadratura de K & las integr de dominio
        # de f. Troceamos con un criterio sencillo
        global nnodos,nelems,nff, nodos,elems,ffx,ffy, ctes,cargas 
        
        trozos= trocear(elems)  # al menos un troceo hacemos
        for i in range (precision_integrar):
            trozos = trocear(trozos)

        # comprobar(trozos) # saca demasiado rollo inutil
        
        if mostrar_integracion:   # info grafica de integracion
            plt.figure('Info: troceo para la integración.')
            
            # dibujo los nodos:
            for i,n in nodos.items():
                x,y = n.x, n.y
                plt.plot(x, y, 'yd', markersize=3)
                a = (np.random.random()-0.5) * anchox/25
                b = (np.random.random()-0.5) * anchoy/25
                plt.text(x+a, y+b, str(i), fontsize=9)

            # dibujo los puntos de Gauss:
            for e in trozos.values():
                for i in range(len(e.x_Gaus[e.tipoe])):
                    psi,eta = e.x_Gaus[e.tipoe][i], e.y_Gaus[e.tipoe][i]
                    x,y = e.dimexy(psi,eta)
                    plt.plot(x,y,'mx', markersize=3)

            # dibujo los bordes en detalle:
            for e in trozos.values():
                e.pinta_elem(colorin='black', ancholin=0.8, tipolin='')

            # ya habra un plt.show() al mostrar resultados

        K, aporta = np.zeros([2*nff,2*nff]), np.zeros([2*nff,1])
        for e in trozos.values():
            for i in range(len(e.x_Gaus[e.tipoe])):
                psi,eta,w = e.x_Gaus[e.tipoe][i], e.y_Gaus[e.tipoe][i], e.w_Gaus[e.tipoe][i]
                x,y = e.dimexy (psi,eta)
                jaco = e.dimeJaco(psi,eta) 
                K += K_integrando (x,y) * w * jaco
                if cargas.volumen:
                    aporta += carga_V_integrando(x,y) * w * jaco
                if cargas.termica:
                    aporta += carga_T_integrando(x,y) * w * jaco

        return(K, aporta)

    def hacer_lambdas():
        # construye como funciones lambda de numpy:
        # la matriz N(xy) de funciones de forma
        # la matriz LN(xy) de derivadas de las ff (se usa en K y en carga T)
        # el integrando de K  (LN)^t * D * (LN)  tambien como func de xy
        # la carga de volumen X[] y la carga termica T (todo de xy)

        # usa las listas ffx, ffy, pero no las modifica (no necesita global)
        # sympy solo se utiliza dentro de  este scope hacer_lambdas

        global nnodos,nelems,nff, nodos,elems,ffx,ffy,N_xy, ctes,cargas 

        x,y = sp.symbols('x,y')
        fforma, lff = sp.zeros(2,2*nff), sp.zeros(3,2*nff)
        
        for i in range (nff):
            # construir la matriz N (osea ff) en entorno sympy
            fforma[0, 2*i]    = sp.sympify(ffx[i])
            fforma[1, 2*i+1]  = sp.sympify(ffy[i])
            
            # construir la matriz LN (osea lff) en entorno sympy
            ffx_x= sp.diff(ffx[i],x)
            ffx_y= sp.diff(ffx[i],y)
            ffy_x= sp.diff(ffy[i],x)
            ffy_y= sp.diff(ffy[i],y)
            lff[0, 2*i]   = ffx_x
            lff[1, 2*i+1] = ffy_y
            lff[2, 2*i]   = ffx_y
            lff[2, 2*i+1] = ffy_x
        
        #print('\nfunciones sympy que entiende. fforma & lff:')
        #print(fforma,'\n', lff)
        
        # construir los integrandos & lambdificarlos

        integrando_K = lff.transpose() * (matrizD * lff)
        K_numeric= sp.lambdify([x,y], integrando_K, 'numpy')
        
        if cargas.volumen:
            Xgrande=sp.zeros(2,1)
            Xgrande[0] = sp.sympify(cargas.volumen[0])
            Xgrande[1] = sp.sympify(cargas.volumen[1])
            integrando_carga_V= fforma.transpose() * Xgrande
            carga_V_numeric = sp.lambdify([x,y], integrando_carga_V, 'numpy')
        else:
            carga_V_numeric = 0.  # que haya algo en el return

        if cargas.termica:
            eps_0 = sp.zeros(3,1)
            if ctes[3] : # tension plana
                eps_0[0] = sp.sympify(cargas.termica[0])*ctes[2]
            else:        # deformacion plana
                eps_0[0] = sp.sympify(cargas.termica[0])*ctes[2]*(1+ctes[1])
            eps_0[1] = eps_0[0]
            eps_0[2] = 0
            integrando_carga_T= lff.transpose() * (matrizD * eps_0)
            carga_T_numeric = sp.lambdify([x,y], integrando_carga_T, 'numpy')
        else:
            carga_T_numeric = 0.

        lff_numeric= sp.lambdify([x,y], lff, 'numpy')
        fforma_numeric= sp.lambdify([x,y], fforma, 'numpy')

        return (fforma_numeric, lff_numeric, K_numeric, carga_V_numeric, carga_T_numeric)

    def apendizar():
        if nfcompleto =='':
            texto = 'Debe guardar el problema antes de usar esta opción.'
            messagebox.showinfo(message=texto,title='MAL',parent=v5)
            return()
        f=open(nfcompleto, 'a')
        print(res_txt, file=f)
        f.close()
        texto ='Los resultados han sido añadidos'
        texto+='al final del fichero de datos:\n'
        texto+= nfcompleto
        messagebox.showinfo(message=texto,title='Hecho',parent=v5)


    comprobar(elems)

    # construir la matriz D
    matrizD = haz_D()
    
    # construir el integrando de la matriz de rigidez (K), de 
    # cargas eps_0 (LN), & para f_pun (N) como lambdas de numpy (!):
    N_xy,LN_xy,K_integrando,carga_V_integrando,carga_T_integrando = hacer_lambdas()
    
    # evaluar las integrales por subdivision y cuadratura
    K_integrada, f_cargas = integrar_en_V()
    
    # tratamiento de las cargas 
    # f_cargas ya viene con las aportaciones de dominio, shape([2*nff,1]), NO [2*nff,] !!
    if cargas.contorno[0]:
        for carga in cargas.contorno[1:]:
            f_cargas += trata_s(carga).reshape(2*nff,1)
    if cargas.puntual[0]:
        for carga in cargas.puntual[1:]:
            f_cargas += trata_p(carga).reshape(2*nff,1)

    # Resuelvo las [a] de la aproximacion (por ahora K sera regular)
    a_despl = np.linalg.solve(K_integrada, f_cargas).reshape(2*nff)
    
    ######################################################################
    # La aproximación esta calculada; las tensiones se calculan al vuelo #
    # Una ventana con texto (y no el terminal) mostrara los resultados.  #
    ######################################################################
    
    res_txt  = '\n'+'#'*75 + '\nIdentificación: {:}\n'.format(nfcompleto)
    res_txt += '\nParametros [a] de la aproximacion:\n'
    res_txt += '  i            a_x                a_y\n'
    
    for i in range(len(a_despl)//2):
        res_txt += '{:4d}   {:+16.6e}   {:+16.6e}\n'.format(
                    i, a_despl[2*i],a_despl[2*i+1])

    # Calculo los desplazamientos nodales
    res_txt += '\nDesplazamientos en los nodos de geometria\n'
    res_txt += ' nodo         ux               uy\n'
    for i,n in nodos.items():
        try:
            a = np.matmul(N_xy(n.x, n.y), a_despl).reshape(2)
        except (ZeroDivisionError):
            res_txt += '\nDivision /0 en el despl. del nodo '+str(i)+ '. Se sustituye por 0.\n'
            a=[0., 0.]
        n.ux, n.uy = a[0], a[1]
        res_txt += '{:3d}  {:+16.6e}  {:+16.6e}\n'.format(i, n.ux, n.uy)

    # Calculo tensiones en los nodos
    if cargas.termica: 
        c_rara=1. if ctes[3] else (1.+ctes[1]) # si es DP es 1+nu en eps0
        eps0_xy= lambda x,y : np.array([eval(cargas.termica[0]), 
                            eval(cargas.termica[0]), 0.]) * ctes[2] * c_rara
    else:
        eps0_xy= lambda x,y : np.zeros(3)

    for i,n in nodos.items():
        try:
            kk= np.matmul(LN_xy(n.x, n.y), a_despl) - eps0_xy(n.x,n.y)
            sigma = np.matmul( matrizD, kk)
        except(ZeroDivisionError):
            res_txt += '\nDivision /0 en la tension del nodo '+str(i)+'. Se sustituye por Inf.\n'
            sigma=[float('Inf'), float('Inf'), float('Inf')]            
        n.sigma = sigma
    # la de von Mises en los nodos la calculo al vuelo de imprimir   
    res_txt += '\nTensiones en los nodos\n'
    res_txt += ' nodo      x          y         sxx          syy          sxy         von_M\n'
    for i,n in nodos.items():
        texto='{:3d}  {:10.5f} {:10.5f}  '.format(i, n.x, n.y)
        for j in range(3): texto += '{:+11.4e}  '.format(n.sigma[j])
        texto += '{:11.4e}'.format(sigma_vm(n.sigma))
        res_txt += texto +'\n'

    # tensiones en polares, si se han solicitado
    if sigma_polar:
        res_txt += '\nTensiones en polares:\n nodo     srr          stt          srt\n'
        for i,n in nodos.items():
            s_pol= a_polares(n.x, n.y, n.sigma[0], n.sigma[1], n.sigma[2])
            res_txt += '{:3d}  {:11.4e}  {:11.4e}  {:11.4e}\n'.format(
                i, s_pol[0], s_pol[1], s_pol[2])
    
    
    v5=Toplevel(v0)
    v5.title('RESULTADOS')
    tcaja = Text(v5, width=80, height=35,wrap='word', font=('Mono',10),
        background='#F0F0F0', foreground='black', border=None, padx=20, pady=12)
    tcaja.grid(column=0, row=0, padx=8, sticky=(N,W,E,S))
    tcaja.insert('1.0',res_txt)

    scb = ttk.Scrollbar(v5, orient=VERTICAL, command=tcaja.yview)
    scb.grid(column=1, row=0, sticky='ns')
    
    tcaja['yscrollcommand'] = scb.set

    apendizar_bt=ttk.Button(v5, text='Guardar con los datos', command=apendizar)
    apendizar_bt.grid(column=0, row=1, pady=4, columnspan=2)
    aviso  = '\nEsto sirve para guardar los resultados anteriores\n'
    aviso += 'en un lugar "lógico" para su consulta posterior:\n'
    aviso += 'se añadirán al fichero de datos, que actualmente es:\n'
    kk= nfcompleto + '\n' if nfcompleto else '      *no hay fichero de datos*\n'
    aviso +=  kk
    aviso += 'ATENCION: si ha modificado el problema y no lo ha\n'
    aviso += 'guardado estará añadiendo resultados que no\n'
    aviso += 'corresponden a los datos.\n'
    
    CreateToolTip(apendizar_bt, aviso)

    tcaja['state']='disabled'

    v5.grid_columnconfigure(0, weight=1)
    v5.grid_rowconfigure(0, weight=4)
    v5.grid_rowconfigure(1, weight=1)
    #v5.geometry('+640+90')
    


    # Resultados graficos
    boton_niveles.state(['!disabled'])
    salida_grafica()
    




#####################################################################
####                    PROGRAMA PRINCIPAL                       ####
#### Basicamente consiste en construir el GUI de control y dejar ####
#### que se llamen desde ahi las ordenes que de el usuario, que  ####
#### se programan como funciones (def) aparte. Lo encuentro      ####
#### menos lioso que programar el gui como función y tener que   ####
#### andar con un monton de globales etc                         ####
#####################################################################

elige, n_f, nfcompleto, res_txt, sin_gui = '', '', '', '', False
#dir_home=path.expanduser('~') # lo hace mejor el solo 
extremosx, extremosy, anchox, anchoy, n_nodos_tipoe = [], [], 0., 0., [8,6]
nnodos,nelems,nff, nodos,elems,ffx,ffy = 0,0,0,{},{},[],[]
a_despl, eps0_xy, sigma_polar = [], [], 0
ctes, cargas, matrizD, N_xy, LN_xy = [], Carga(), [],  0., 0.
hard_nnodos, hard_nelems, hard_nff = 24, 6, 24
mostrar_integracion, precision_integrar, figuras_bien, v5 = 0, 0, 0 ,None
niveles_sigma=[['',''],['',''],['',''],['',''],['',''],['',''],['','']]
niveles_sigma_bak=[['',''],['',''],['',''],['',''],['',''],['',''],['','']]


v0=Tk()
presenta_elige() # v0 es auto destruida 

v0=Tk()
v0.title("GalerkinAppr v0.7 - ejemplo")



estilo = ttk.Style()
estilo.configure('jc.TButton', foreground='green')
estilo.configure('jc.TLabelframe.Label', foreground ='green')
estilo.configure('jc_blue.TButton',background='#CADDF3')# ,foreground='#2A6AE6'
estilo.configure('jc_red.TButton',foreground='#9B2803')



frame_general = ttk.Frame(v0, padding='4')
frame_general.grid(sticky=(N, W, E, S))

frame_acciones= ttk.Labelframe(frame_general, text='Acciones',
    style='jc.TLabelframe', width=170,height=550)
frame_acciones.grid(column=0,row=0,ipadx=3)

frame_pestanas= ttk.Notebook(frame_general)
frame_pestanas.grid(column=1, row=0)

frame_dominio=    ttk.Frame(frame_pestanas)
frame_funciones=ttk.Frame(frame_pestanas)
frame_cargas=   ttk.Frame(frame_pestanas)

frame_pestanas.add(frame_dominio, text='Dominio')
frame_pestanas.add(frame_funciones, text='Funciones')
frame_pestanas.add(frame_cargas, text='Cargas')


# botones de acciones
boton_cargar=ttk.Button(frame_acciones, text='cargar', command= a_cargar)
boton_cargar.grid(column=0, row=0, padx=3, pady=(20,5))

boton_calcula=ttk.Button(frame_acciones,text='calcular',
    style='jc_blue.TButton', command= calcula)
boton_calcula.grid(column=0, row=1, padx=3, pady=5) 

boton_guardar=ttk.Button(frame_acciones, text='guardar', command=a_guardar)
boton_guardar.grid(column=0, row=2, padx=3, pady=5)

ttk.Separator(frame_acciones,orient=HORIZONTAL).grid(
    row=3,column=0,pady=20, sticky='ew')

boton_salir=ttk.Button(frame_acciones, text='-salir-', style='jc_red.TButton',
    command=a_salir)
boton_salir.grid(column=0, row=4, padx=3, pady=5)

boton_obviar=ttk.Button(frame_acciones, text='obviar GUI', 
style='jc_red.TButton', command=obviar_gui)
boton_obviar.grid(column=0,row=5, padx=3, pady=5)

boton_limpiar=ttk.Button(frame_acciones,text=' limpiar', style='jc_red.TButton',
    command= a_limpiar)
boton_limpiar.grid(column=0, row=6, padx=3, pady=5)

ttk.Separator(frame_acciones,orient=HORIZONTAL).grid(
    row=7,column=0,pady=20, sticky='ew')

boton_flip_f=ttk.Button(frame_acciones, text='fx<>fy', command= flip_f)
boton_flip_f.grid(column=0, row=8, padx=3, pady=5)
texto = 'Intercambia las funciones de ux con las de uy.\n'
texto +='Esto permite usar las funciones de un problema\n' 
texto +='simétrico para uno antisimétrico o viceversa.\n' 
CreateToolTip(boton_flip_f, texto)

boton_pintaff=ttk.Button(frame_acciones, text='ver ff', command=ver_ff )
boton_pintaff.grid(column=0, row=9, padx=3, pady=5)

ttk.Separator(frame_acciones,orient=HORIZONTAL).grid(
    row=10,column=0, pady=20, sticky='ew')

boton_niveles=ttk.Button(frame_acciones,text='niveles',
    command= niveles_color)
boton_niveles.grid(column=0, row=11, padx=3, pady=5)
texto = 'Permite ajustar el rango de tensiones\n'
texto +='cubierto por la progresión de colores\n'
texto +='en las gráficas (tras "calcular").\n'
CreateToolTip(boton_niveles, texto)
boton_niveles.state(['disabled'])


boton_cerrargr=ttk.Button(frame_acciones,text=' cerrar\ngráficos',
    command= cierraplots)
boton_cerrargr.grid(column=0, row=12, padx=3, pady=5)

boton_opciones=ttk.Button(frame_acciones, text='opciones', command=entraOpc)
boton_opciones.grid(column=0, row=13, padx=3, pady=(5,20))

############################################
# rellenar la ventana-pestana frame_dominio
############################################

frame_nodos=ttk.Labelframe(frame_dominio, text='Nodos', style='jc.TLabelframe')
frame_elems=ttk.Labelframe(frame_dominio, text='Elementos', style='jc.TLabelframe')
frame_ctes= ttk.Labelframe(frame_dominio, text='Material', style='jc.TLabelframe')

frame_nodos.grid(row=0, column=0, rowspan=2, padx=16)
frame_elems.grid(row=0, column=1, sticky='n')
frame_ctes.grid(row=1,column=1, sticky='n')


    # subventana frame_dominio -> frame_nodos 

label_nodos1=ttk.Label(frame_nodos,text='nº').grid(column=0, row=0, padx=3, pady=2)
label_nodos2=ttk.Label(frame_nodos,text='x').grid(column=1, row=0, padx=3, pady=2)
label_nodos3=ttk.Label(frame_nodos,text='y').grid(column=2, row=0, padx=3, pady=2)

filas_nodos = []
for ifilas in range(hard_nnodos):
    fila=[]
    for j in range(3):
        (ancho,a) = (3,4) if j==0 else (14,0)
        e=ttk.Entry(frame_nodos, width=ancho)
        e.grid(row=ifilas+1, column=j, padx=a)
        fila.append(e)
    filas_nodos.append(fila)

    # subventana frame_dominio -> frame_elems
        
label_el1=ttk.Label(frame_elems,text='nº').grid(column=0, row=0, padx=3, pady=2)
label_el2=ttk.Label(frame_elems,text='nodos del elemento').grid(
    column=1, row=0, padx=3, pady=2, columnspan=8)

filas_elems = []
for ifilas in range(hard_nelems):
    fila=[]
    for j in range(9):
        (ancho,a) = (3,12) if j==0 else (5,1)
        e=ttk.Entry(frame_elems, width=ancho)
        e.grid(row=ifilas+1, column=j, padx=(1,a) )
        fila.append(e)
    filas_elems.append(fila)

boton_dibuja =ttk.Button(frame_elems,text='Visualizar', command=pinta_comprobar)
boton_dibuja.grid(column=1, row=hard_nelems+2, columnspan=8, pady=10)

    # subventana frame_dominio -> frame_ctes

label_young=ttk.Label(frame_ctes,text='Módulo de Young').grid(
    column=0, row=0, padx=3, pady=2, sticky='e')
label_poiss=ttk.Label(frame_ctes,text='Coef. Poisson').grid(
    column=0, row=1, padx=3, pady=2, sticky='e')
label_dilat=ttk.Label(frame_ctes,text='Coef. dilatación térmica').grid(
    column=0, row=2, padx=3, pady=2, sticky='e')
label_tp=ttk.Label(frame_ctes,text='Tensión plana').grid(
    column=0, row=3, padx=3, pady=2, sticky='e')
label_dp=ttk.Label(frame_ctes,text='Deformación plana').grid(
    column=0, row=4, padx=3, pady=2, sticky='e')

entry_young=ttk.Entry(frame_ctes, width=10)
entry_young.grid(column=1, row=0, padx=3, pady=2, sticky='w')
entry_poiss=ttk.Entry(frame_ctes, width=10)
entry_poiss.grid(column=1, row=1, padx=3, pady=2, sticky='w')
entry_dilat=ttk.Entry(frame_ctes, width=10)
entry_dilat.grid(column=1, row=2, padx=3, pady=2, sticky='w')

tens_plana=StringVar()
tens_plana.set('1')
t_plana=Radiobutton(frame_ctes, variable=tens_plana, value='1',background='#D9D9D9') 
d_plana=Radiobutton(frame_ctes, variable=tens_plana, value='0',background='#D9D9D9')  #DCE8D4
t_plana.grid(column=1, row=3)
d_plana.grid(column=1, row=4)



#################################################
# rellenar la ventana-pestana de frame_funciones
#################################################

label_n=ttk.Label(frame_funciones,text='nº').grid(
    column=0, row=0, padx=3, pady=2)
label_paux=ttk.Label(frame_funciones,text='funciones para ux').grid(
    column=1, row=0, padx=3, pady=2)
label_pauy=ttk.Label(frame_funciones,text='funciones para uy').grid(
    column=2, row=0, padx=3, pady=2)

filas_funcs=[]
for ifilas in range(hard_nff):
    fila=[]
    for j in range(3):
        ancho=3 if j==0 else 40
        e=ttk.Entry(frame_funciones, width=ancho)
        e.grid(row=ifilas+1, column=j, padx=3)
        if j==0:
            e.insert(0, str(ifilas))
            e.state(['disabled'])
        fila.append(e)
    filas_funcs.append(fila)


#################################################
# rellenar la ventana-pestana de frame_cargas
#################################################

frame_xraya=ttk.Labelframe(frame_cargas, text='Cargas distribuidas de contorno', 
            style='jc.TLabelframe')
frame_fuerzas=ttk.Labelframe(frame_cargas, text='Cargas puntuales', style='jc.TLabelframe')
frame_xgrande= ttk.Labelframe(frame_cargas, text='Carga de volumen', style='jc.TLabelframe')
frame_temperatura= ttk.Labelframe(frame_cargas, text='Carga térmica', style='jc.TLabelframe')

frame_xraya.grid(row=0, column=0, padx=(10,0), pady=(16,6), sticky='w')
frame_fuerzas.grid(row=1, column=0, padx=(10,0), pady=10, sticky='w')
frame_xgrande.grid(row=2,column=0, padx=(10,0), pady=10, sticky='w')
frame_temperatura.grid(row=3,column=0, padx=(10,0), pady=10, sticky='w')

    # subventana frame_cargas -> frame_xraya

kk=ttk.Label(frame_xraya, text='¿coor nt?').grid(row=0,column=0)
kk=ttk.Label(frame_xraya, text='nodos').grid(row=0,column=1, columnspan=3)
kk=ttk.Label(frame_xraya, text='componentes_1').grid(row=0,column=4,columnspan=3)
kk=ttk.Label(frame_xraya, text='componentes_2').grid(row=0,column=7,columnspan=3)

filas_xraya=[]
espacios=[ 0,0, (0,8), 0,0, (0,8), 0,0, (0,6)]
for ifilas in range(6): # sitio para 6 cargas de contorno 
    fila=[1]
    fila[0]=StringVar()
    e=ttk.Checkbutton(frame_xraya, variable=fila[0])
    e.grid(column=0, row=ifilas+1)
    for j in range(9):
        ancho= 3 if j<3 else 8
        e=ttk.Entry(frame_xraya, width=ancho)
        e.grid(row=ifilas+1, column=j+1, padx=espacios[j])
        fila.append(e)
    filas_xraya.append(fila)
    
    # subventana frame_cargas -> frame_fuerzas

kk=ttk.Label(frame_fuerzas, text='nodo').grid(row=0, column=0, sticky='e')
kk=ttk.Label(frame_fuerzas, text='componente x').grid(row=1, column=0)
kk=ttk.Label(frame_fuerzas, text='componente y').grid(row=2, column=0)

filas_fuerzas=[]
for ifilas in range(6): # en realidad lo dispongo en columnas // sitio para 6fuerzas
    fila=[]
    for j in range (3):
        ancho= 3 if j==0 else 8
        e=ttk.Entry(frame_fuerzas, width=ancho)
        e.grid(row=j, column=ifilas+1, padx=12)
        fila.append(e)
    filas_fuerzas.append(fila) # insisto, en el gui esta en columnas


    # subventana frame_cargas -> frame_xgrande

kk=ttk.Label(frame_xgrande, text='componente x: ').grid(row=0, column=0)
kk=ttk.Label(frame_xgrande, text='componente y: ').grid(row=1, column=0)

filas_xgrande=[]
for i in range(2):
    e=ttk.Entry(frame_xgrande, width=60)
    e.grid(row=i, column=1, padx=3)
    filas_xgrande.append(e)
    
    # subventana frame_cargas -> frame_temperatura

kk=ttk.Label(frame_temperatura, text='Temperatura(x,y): ').grid(row=0, column=0)

filas_temperatura=[]
e=ttk.Entry(frame_temperatura, width=57)
e.grid(row=0, column=1, padx=3)
filas_temperatura.append(e)


# el GUI esta. Si se pidio el prb de ejemplo, rellenar con sus datos y hacer
if elige=='default': 
    rellena_default()
    rollete_default()
    calcula()

v0.protocol('WM_DELETE_WINDOW', a_salir)
v0.mainloop()




