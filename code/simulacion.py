import numpy as np
import matplotlib.pyplot as plt
import ROOT as root
import scipy as sp

# Funcion Uniformidad Integral.
# Input: N_pix, C_med.
# Output: unif_integral (escalar)



def field_of_view(N_pix, filas_eliminadas, columnas_eliminadas):
    fov = [N_pix - filas_eliminadas, N_pix - columnas_eliminadas]
    return fov


def matriz_uniform(N_pix, C_med, field_of_view):                                                   
    matriz_cuentas = np.zeros((N_pix,N_pix))
    filas_ini = int(N_pix - field_of_view[0])
    filas_fin = int(field_of_view[0])
    col_ini = int(N_pix - field_of_view[1])
    col_fin =int(field_of_view[1])
    
    for i in range(filas_ini, filas_fin):
        for j in range(col_ini, col_fin):
            matriz_cuentas[i,j] = np.random.normal(C_med, np.sqrt(C_med))
    return matriz_cuentas


def filtro_9p(matriz_cuentas,N_pix,filtro=1):
    matriz_filtrada = np.zeros((N_pix,N_pix))
    
    for i in range(4,N_pix-4):
        for j in range(12,N_pix-12):
            P1=matriz_cuentas[i-1,j-1]
            P2=matriz_cuentas[i-1,j]
            P3=matriz_cuentas[i-1,j+1]
            P4=matriz_cuentas[i,j-1]
            P5=matriz_cuentas[i,j]
            P6=matriz_cuentas[i,j+1]
            P7=matriz_cuentas[i+1,j-1]
            P8=matriz_cuentas[i+1,j]
            P9=matriz_cuentas[i+1,j+1]
            if(filtro==1):
                matriz_filtrada[i,j]=(P1*1 + P2*2 + P3*1 +2*P4 + 4*P5 + 2*P6 + 1*P7 + 2*P8 + 1*P9)/(16)
            if(filtro==0):
                matriz_filtrada[i,j]=P5
    
    "Aqui extraigo la submatriz para continuar con el analisis"
    start_row = 4
    size_row = N_pix-8
    start_col=12
    size_col=N_pix-24
    submatrix = matriz_filtrada[start_row:start_row + size_row, start_col:start_col + size_col]
    return submatrix


    


def uniformidad_integral(matriz_cuentas):
    max = np.max(matriz_cuentas)
    min = np.min(matriz_cuentas)
    unif_integral = 100*(max - min)/(max + min)
    return unif_integral

def analisis(N_pix, C_med,filas_eliminadas, columnas_eliminadas,filtro=1):
    fov = field_of_view(N_pix, filas_eliminadas, columnas_eliminadas)
    matriz_cuentas = matriz_uniform(N_pix,C_med,fov)
    matriz_filtrada = filtro_9p(matriz_cuentas,N_pix,filtro)
    unif_integral = uniformidad_integral(matriz_filtrada)
    unif_integral = uniformidad_integral(matriz_cuentas)
    return unif_integral
    


#Generar distribucion Histograma sin defecto
def dist(N_pix,C_med,num_cuentas,filas_eliminadas,columnas_eliminadas,filtro=1):
    d = np.zeros(num_cuentas)
    for i in range(num_cuentas):
        d[i] = analisis(N_pix, C_med,filas_eliminadas, columnas_eliminadas,filtro)
    distribucion = np.sort(d)
    return distribucion



#Introduciendo los defectos en las matrices sin defectos 
def matriz_defectos(t,k,N_pix,C_med,fov):
    matriz = matriz_uniform(N_pix,C_med,fov)
    n=int(t*(N_pix/64))
    for i in range(int(N_pix/2),int(N_pix/2)+n):
        for j in range(int(N_pix/2),int(N_pix/2)+n):
             matriz[i,j] = np.random.normal((1+k/100)*C_med, np.sqrt((1+k/100)*C_med))

    return matriz

def analisis_defectos(t,k,N_pix, C_med,filas_eliminadas, columnas_eliminadas,filtro=1):
    fov = field_of_view(N_pix, filas_eliminadas, columnas_eliminadas)
    matriz_cuentas = matriz_defectos(t,k,N_pix,C_med,fov)
    matriz_filtrada = filtro_9p(matriz_cuentas,N_pix,filtro)
    unif_integral = uniformidad_integral(matriz_filtrada)
    return unif_integral

def dist_defectos(t,k,N_pix,C_med,num_cuentas,filas_eliminadas,columnas_eliminadas,filtro=1):
    d = np.zeros(num_cuentas)
    for i in range(num_cuentas):
        d[i] = analisis_defectos(t,k,N_pix,C_med,filas_eliminadas,columnas_eliminadas,filtro)
    distribucion = np.sort(d)
    return distribucion

#incertidumbres de un histograma
"400 bins a 0.025 de 0.0 a 10.0"

def incertidumbre(histograma,num_cuentas): #histograma debe estar normalizado
    freq=np.zeros(400)
    sigma=np.zeros(400)
    for i in range(400):
        freq[i]=histograma.GetBinContent(i)
        if freq[i]<1.0:
           sigma[i]=np.sqrt((1/num_cuentas)*freq[i]*(1-freq[i]))
        else:
            sigma[i]=0.0
    return sigma

def name(t,k,N_pix,nct,C_med):
    return "d_"+"t="+str(t)+"_k="+str(k)+"_"+str(N_pix)+"P_"+str(nct)+"Mc_"+"Cmed_"+str(C_med)

def in_defectos(t,k,C_med,matriz_sana):
    N_pix = np.shape(matriz_sana)[0]
    matriz = matriz_sana
    n=int(t*(N_pix/64))
    for i in range(int(N_pix/2),int(N_pix/2)+n):
        for j in range(int(N_pix/2),int(N_pix/2)+n):
             matriz[i,j] = np.random.normal((1+k/100)*C_med, np.sqrt((1+k/100)*C_med))
    return matriz

'''funcion para nombrar a los hitogramas: xxx(num_pixel)-xx(mill.cuentas)-fx(filtro, 0=no filtro)-xx(tamaño, 0=no hay defecto)-xx(contraste*10)
    nombre = 064-80-f1-2-10 significa 64pix 80Mc filtro9p t=2 k=1.0
    Cmed no es un parametro independiente, se calcula a partir de 80Mc y 64pix'''

def histogram_name(N_pix,nct,filtro=1,t=0,k=0):
    N_pix = str(int(N_pix))
    if len(N_pix)<3:
        N_pix = '0'+N_pix
    
    nct = str(int(nct))
    if (nct=="5"): nct="0"+nct
    filtro = str(int(filtro))
    
    if t==0:k=0
    if k==0:t=0
    t = str(t)
    if int(k)==10:
        k="100"
    else:
        k= str(int(k*10))
    name = N_pix +'-'+nct+ '-f' + filtro + '-' + t + '-' + k + ".root"
    return name


'''
Funciones para el ANALISIS ESTADISTICO, incluyendo::
-- Integracion
-- Generar curvas ROC, dibujar curvas ROC
-- Generar curvas area contraste, dibujar curvas area contraste
-- Generar curvas contraste detalle, dibujar curvas contraste detalle
'''

#INTEGRAR 




def trapecio(x,y):
    integral = 0
    for i in range(len(x)-1):
        delta_x = x[i+1]-x[i]
        integral = integral + 0.5*(y[i+1]+y[i])*delta_x
    return abs(integral)

def sorteo_montecarlo(x,y,xe,ye):
    #Inicializo los nuevos puntos x,y a sortear
    x_mc = np.zeros(len(x))
    y_mc = np.zeros(len(x))
    for i in range(len(x)):
        x_mc[i]=np.random.normal(x[i],xe[i])
        y_mc[i]=np.random.normal(y[i],ye[i])
    #Uso x_errors, y_errors inicializados a cero porque no influyen
    return x_mc,y_mc




def montecarlo_trapecio(x,y,xe,ye,error_lineal=False):
    integral_inicial = trapecio(x,y)
    integral = np.array([])

    if(error_lineal==False):
            for i in range(1000):
                x_mc,y_mc=sorteo_montecarlo(x,y,xe,ye)
                integral = np.append(integral,trapecio(x_mc,y_mc)) #Da igual que valores de xe,ye tomememos. De hecho creo que lo voy a borrar
            integral = np.sort(integral)
            sigma = np.std(integral)
    if(error_lineal==True):
        sigma2=0
        for i in range(len(x)-1):
            delta_x = x[i+1]-x[i]
            delta_y = 0.5*(y[i+1]-y[i])
            sigma_x = np.sqrt(xe[i]**2+xe[i+1]**2)
            sigma_y = 0.5*np.sqrt(ye[i]**2+ye[i+1]**2)
            sigma2 = sigma2 + (delta_y*sigma_x)**2 + (delta_x*sigma_y)**2
        sigma = np.sqrt(sigma2)
    return integral_inicial,sigma


#--------------------CURVAS ROC y AREA CONTRASTE------------------------
''' Voy a crear dos funciones principales:
1. curvaroc(n_pix,nct,t,k)
Input: modo de adquisicion
Output: Datos(x,y,xe,ye) y gráfica curva roc

2. areacontraste(n_pix,nct,t)
Input: modo de aquisicion salvo k (ahora es una variable)
Output: Objeto TGraphErrors y gráfica curva area contraste
---------de momento no realizo en ajuste dentro de esta funcion.

Tambien he incluido las siguientes funciones secundarias:
3. ajuste_ac(x,y,xe,ye)
Input: valores de la curva area-contraste
Output: parametros a,b del ajuste a la curva logistica
4. dibujarroc(x,y,xe,ye)
5. dibujarAC(x,y,xe,ye) ---- sirven para dibujar y guardar en un archivo las curvas roc y area-contraste
'''

def curvaroc(N_pix,nct,t,k,dibujar=False,filtro=1): #de momento el filtro es f1

    file_path = "./histogramas/64x64/"
    name_sana = histogram_name(N_pix,nct,filtro)
    name_defec = histogram_name(N_pix,nct,filtro,t,k)
    file_sana = root.TFile(file_path+name_sana)
    file_defec = root.TFile(file_path+name_defec)
    hist_sana = file_sana.Get(name_sana)
    hist_defec = file_defec.Get(name_defec)
    sana = hist_sana.GetCumulative()
    defec = hist_defec.GetCumulative()
    nbins = hist_sana.GetNbinsX() #nbins==400 Se debería aumentar?
    x,y = np.zeros(400),np.zeros(400)
    xe,ye = incertidumbre(sana,10000),incertidumbre(defec,10000)
    for i in range(400):
        y[i] = 1-defec.GetBinContent(i)
        x[i] = 1-sana.GetBinContent(i)
    if dibujar==True: dibujarroc(x,y,xe,ye,root.kBlue)
    #se podria plantear la opcion de devolver un .cv con los datos
    return x,y,xe,ye

def areacontraste(N_pix,nct,t,dibujar=False,lineal_error=False): #k no es input, es una variable de output
    k_list=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
    integral_list=[]
    error_list=[]
    for k in k_list:
        x,y,xe,ye = curvaroc(N_pix,nct,t,k)
        integral,error = montecarlo_trapecio(x,y,xe,ye,lineal_error)
        integral_list.append(integral)
        error_list.append(error)
    x = np.array(k_list)
    xe = np.zeros(len(x))
    y = np.array(integral_list)
    ye = np.array(error_list)
    #dibujarAC(x,y,xe,ye,root.kRed)
    return x,y,xe,ye
    
def ajuste_ac(x,y):
    #0. Variable auxiliar para decidir si el ajuste se ha podido realizar
    fitted = False
    #1. Transformar los datos para realizar un ajuste lineal.
    x_fit,y_fit = np.zeros(len(x)),np.zeros(len(x))
    for i in range(len(x)):
        if 0.51<y[i]<0.99: #No consideros los puntos en los extremos saturados
            fitted = True
            x_fit[i] = x[i] #x no cambia
            y_fit[i] = np.log(1/(2*y[i]-1)-1) #funcion de ajuste invertida -- falta considerar la transformacion de los errores
    if fitted == True:
        x_fit = x_fit[x_fit!=0] #Elimino los elementos no nulos
        y_fit = y_fit[y_fit!=0]
        x_min = np.min(x_fit)
        x_max = np.max(x_fit)
        curvalineal = root.TGraph(len(x_fit),x_fit,y_fit)
        linearfit = root.TF1("linearfit","1++x",x_min,x_max) # a+b*x
        curvalineal.Fit(linearfit,"C")
        linearfit = curvalineal.GetFunction("linearfit")
        a = linearfit.GetParameter(0)
        b = linearfit.GetParameter(1)
        return a,b, fitted
    if fitted == False:
        x0 = 0
        m = 0
        for i in range(len(x)-1):
            if y[i+1]-y[i]>0.2:
                x0 = x[i]
                m = (y[i+1]-y[i])/(x[i+1]-x[i])
            return x0, m, fitted
    

    
    
  
def dibujarroc(x,y,xe,ye,color): #los argumentos son los valores de fvn y ffn; color es un elemento kcolor
    roc = root.TGraphErrors(400)
    for i in range(400):
        roc.SetPoint(i,x[i],y[i])
        roc.SetPointError(i,xe[i],ye[i])
    roc.SetTitle("Receiver Operating Characteristic (ROC)")
    roc.GetXaxis().SetTitle("False Positive Rate")
    roc.GetYaxis().SetTitle("True Positive Rate")
    canvas = root.TCanvas('','',800,600)   
    roc.SetMarkerColor(color)
    roc.SetLineColor(color)
    roc.SetMarkerStyle(24)
    roc.Draw("APX") #los errores de momento son menores al tamaño de Marker
    canvas.SaveAs("roc_plot.png")

def dibujarAC(x,y,xe,ye,color): #los argumentos son los valores de fvn y ffn; color es un elemento kcolor
    roc = root.TGraphErrors(len(x))
    for i in range(len(x)):
        roc.SetPoint(i,x[i],y[i])
        roc.SetPointError(i,xe[i],ye[i])
    ajuste = root.TGraph(10000)
    x_fit = np.zeros(10000)
    y_fit = np.zeros(10000)
    a,b,fitted = ajuste_ac(x,y)
    for i in range(10000):
        x_fit[i] = i*10/10000
        y_fit[i] = 0.5*(1+1/(1+np.exp(a+b*x_fit[i])))
        ajuste.SetPoint(i,x_fit[i],y_fit[i])
    ajuste.SetLineColor(color-5)
    roc.SetTitle("Curva AC: Area-Contraste")
    roc.GetXaxis().SetTitle("k - contraste")
    roc.GetYaxis().SetTitle("Aroc - Area curva roc")
    canvas = root.TCanvas('','',800,600)   
    roc.SetMarkerColor(color)
    roc.SetLineColor(color)
    roc.SetMarkerStyle(24)
    roc.Draw("AP") #los errores de momento son menores al tamaño del Marker
    ajuste.Draw("SAME")
    canvas.SaveAs("ac_plot.png")

#-------------------CURVAS CONTRASTE DETALLE-------------------
def contrastedetalle(Npix,nct):
    t_list = np.array([1,2,3,4,5])
    k_list = np.zeros(len(t_list))
    sigma_list = np.zeros(len(t_list))
    for i in range(len(t_list)):
        x,y,xe,ye = areacontraste(Npix,nct,t_list[i],False)
        a,b,fitted=ajuste_ac(x,y)
        if fitted == True:
            aroc = 0.80 #Esto puede ser modificado a 0.8 -- 0.9
            if -0.01<b<0.01: break
            #k = (1/b) * (-a + np.log(1/(2*aroc-1)-1))
            k = (1/b) * (-a + np.log(1/(2*aroc-1)-1))
            #implemento el montecarlo directamente
            montecarlo=np.zeros(2000)
    
            for j in range(2000):
                ymc = np.zeros(len(x))
                for l in range(len(x)):
                    ymc[l]=np.random.normal(y[l],ye[l])
                amc,bmc,fitted=ajuste_ac(x,ymc)
                if abs(bmc)>0.01:
                    kmc = (1/bmc)*(-amc + np.log(1/(2*aroc-1)-1))
                montecarlo[j] = kmc
            k_list[i] = k
            sigma_list[i] = np.std(montecarlo)
        if fitted == False:
            x0 = a
            m = b
            k = (0.8-0.5)/m + x0
            k_list[i] = k
            sigma_list[i] = 0.0

    return t_list,k_list,sigma_list


def dibujarCD(x,y,xe,ye,color,nct): #los argumentos son los valores de fvn y ffn; color es un elemento kcolor
    roc = root.TGraphErrors(len(x))
    for i in range(len(x)):
        roc.SetPoint(i,x[i],y[i])
        roc.SetPointError(i,xe[i],ye[i])
    roc.SetTitle("Curva CD: Contraste-Detalle")
    roc.GetXaxis().SetTitle("t^{2} - tamano del defecto (detalle)")
    roc.GetYaxis().SetTitle("k_50 - contraste")
    canvas = root.TCanvas('','',800,600)   
    roc.SetMarkerColor(color)
    roc.SetLineColor(color)
    roc.SetMarkerStyle(24)
    roc.Draw("AP") #los errores de momento son menores al tamaño del Marker
    canvas.SaveAs("cd-"+str(nct)+".png")