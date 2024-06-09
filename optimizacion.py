# Francisco Sciaroni
import matplotlib.pyplot as plt
import numpy as np
import random

class NewtonMethod():
    def find_root(self, func, derivative, x0 = 0.01, x_tolerance=1e-5, y_tolerance = 1e-5, max_iterations=100000, x_step = 1_000):
        x = x0
        for i in range(max_iterations):
            f_x = func(x)
            df_x = derivative(x)
            
            if abs(f_x) < y_tolerance:
                return x

            if df_x == 0:
                raise ValueError(f"Derivative is zero at x = {x}, cannot continue Newton's method")
            
            x_new = x - f_x / df_x
            if abs(x_new - x) < x_tolerance:
                return self.find_root(func, derivative, x0 = random.uniform(x0-x_step, x0+x_step), 
                                      x_tolerance = x_tolerance, y_tolerance = y_tolerance, 
                                      max_iterations = max_iterations - i - 1, x_step = x_step)
            x = x_new

        return None

class Myarray():
    def __init__(self,lista,r,c, by_row = True):
        self.rows = r
        self.columns = c
        self.elems = lista
        self.orden = by_row
        if self.orden == False:
            self.orden = True
            self.elems = self.Orden().elems            
    
    def Orden(self): 
        new_elems = []
        for i in range(self.rows):
            row = self.elems[i::self.rows]
            for elem in row:
                new_elems.append(elem)
        return Myarray(new_elems, self.rows, self.columns, True)
    
    def __str__(self):
        matrix_str = ""
        for i in range(self.rows):
            for j in range(self.columns):
                matrix_str += str(self.GetElem(i+1, j+1)) + " "  
            matrix_str += "\n"
        return matrix_str
    
    def Identidad(self,j,k):
        if j != k:
            raise ValueError('Identity must be a square matrix')
        else:
            elems = [0]*(j*k)
            salida = Myarray(elems, j, k)        
            for i in range(len(elems)):
                coord = salida.GetCoords(i)
                if coord[0] == coord[1]:
                    elems[i] = 1                    
            return salida
    
    def __add__(self, b): 
        if isinstance(b, int) or isinstance(b,float):
            self.elems = [elem + b for elem in self.elems]
            salida = Myarray(self.elems, self.rows, self.columns)
        
        else:
            if not (self.rows == b.rows) and (self.columns == b.columns):
                raise ValueError('the number of rows must be equal to the number of columns')
            new_elems = []
            for i in range(len(self.elems)):
                new_elems.append(self.elems[i]+b.elems[i])
            salida = Myarray(new_elems, self.rows, self.columns)
        return salida
    
    def __radd__(self, b):
        return self.__add__(b)
    
    def __sub__(self,b): 
        if isinstance(b, int) or isinstance(b,float):
            self.elems = [elem - b for elem in self.elems]
            salida = Myarray(self.elems, self.rows, self.columns)
        else:
            if not (self.rows == b.rows) and (self.columns == b.columns):
                raise ValueError('the number of rows must be equal to the number of columns')
            new_elems = []
            for i in range(len(self.elems)):
                new_elems.append(self.elems[i]-b.elems[i])
            salida = Myarray(new_elems, self.rows, self.columns)        
        return salida
    
    def __rsub__(self, b):
        return self.__sub__(b)
    
    def __mul__(self, mult): #arreglar para cuando la otra matriz False
        new_elems = []
        if isinstance(mult, int) or isinstance(mult,float):
            for i in range(len(self.elems)):
                new_elems.append(self.elems[i]*mult)            
            salida = Myarray(new_elems, self.rows, self.columns)
            
        else:
            if self.columns != mult.rows:
                raise ValueError('the number of rows of array A must be equal to the number of columns in array B')
            for i in range(self.rows):
                fila = self.GetRow(i+1)
                for w in range(mult.columns):
                    col = mult.GetCol(w+1)
                    suma = 0
                    for x in range(len(fila)):
                        suma += (fila[x]*col[x])
                    new_elems.append(suma)                    
            salida = Myarray(new_elems, self.rows, mult.columns) 
        return salida
    
    def __rmul__(self, mult):
        return self.__mul__(mult)
    
    def __pow__(self,b):
        if isinstance(b, int) or isinstance(b,float):
            self.elems = [elem **b for elem in self.elems]
            salida = Myarray(self.elems, self.rows, self.columns)
        else:
            raise ValueError('the pow must be an int or float')    
        return salida
    
    def GetPos(self, j, k):
        check = j <= self.rows and k <= self.columns        
        if check == False:
            raise ValueError('rows or columns out of range')
        
        else:
            div = int(len(self.elems)/self.rows)
            index = div*(j-1)+k-1           
            return index 
        
       
    def GetCoords(self, m):
        check = m <= len(self.elems)-1               
        if check == False:
            raise ValueError('index out of range')
        
        else:
            ubicacion = m + 1 
            div = int(len(self.elems)/self.rows)                                
            cont = 0
            row = 1
            col = 1
            for i in range(1, len(self.elems)+1):       
                if ubicacion > i:
                    cont += 1
                    col += 1
                    if cont == div:
                        row += 1
                        col = 1
                        cont = 0
            salida = (row, col)
            
            return salida
                            
    def Switch(self):  
        new_elems = []
        if self.orden == True:
            for k in range(1,self.columns + 1):
                for j in range(1,self.rows + 1):
                    new_elems.append(self.elems[self.GetPos(j,k)])
        else:
            for j in range(1,self.rows + 1):
                for k in range(1,self.columns + 1):
                    new_elems.append(self.elems[self.GetPos(j,k)])  
                        
        return Myarray(self.elems, self.rows, self.columns, not (self.orden))
        
   
    def GetRow(self, j): 
        check = j <= self.rows              
        if check == False:
            raise ValueError('index out of range')
        
        else:
            div = int(len(self.elems)/self.rows)
            row = self.elems[div*(j-1):div*j]
            return row         
        
    def GetCol(self, k):
        check = k <= self.columns              
        if check == False:
            raise ValueError('index out of range')
        
        else:
            column = self.elems[k-1::self.columns]
            return column
            
        
    def GetElem(self,j,k):
        check = j <= self.rows and k <= self.columns        
        if check == False:
            raise ValueError('rows or columns out of range')
        
        else:
            indice = self.GetPos(j,k)
            elemento = self.elems[indice]                
            return elemento
            

    def GetSubmatrix(self, row_list, col_list):
        if not (isinstance(row_list, list) and isinstance(col_list, list)):
            raise ValueError('rowslist and colslist should be lists')
        
        new_mat = self
        cant_rows_eliminadas = 0
        for row in row_list:
            new_mat = new_mat.DelRow(row - cant_rows_eliminadas)
            cant_rows_eliminadas += 1
        
        cant_cols_eliminadas = 0
        for col in col_list:
            new_mat = new_mat.DelCol(col - cant_cols_eliminadas)
            cant_cols_eliminadas += 1
                    
        return new_mat
        

    def DelRow(self,j):
        check = self.rows >= j
        if check == False:
            raise ValueError ('index out of range')
        else:
            new_mat = []        
            for i in range(len(self.elems)):
                coord = self.GetCoords(i)
                x_coord = coord[0]
                if x_coord != j:
                    new_mat.append(self.elems[i])
        
            return Myarray(new_mat, self.rows-1, self.columns,self.orden)
    
    def DelRow2(self,j):
        check = self.columns >= j
        if check == False:
            raise ValueError ('index out of range')
        else:
            iden = self.Identidad(self.rows, self.columns)
            iden = iden.DelRow(j)            
            return iden*self
    
        
    def DelCol(self,k):
        check = self.columns >= k
        if check == False:
            raise ValueError ('index out of range')
        else:
            new_mat = []
            for i in range(len(self.elems)):
                coord = self.GetCoords(i)
                y_coord = coord[1]
                if y_coord != k:
                    new_mat.append(self.elems[i])
            
            return Myarray(new_mat,self.rows,self.columns-1,self.orden)
    
    def DelCol2(self,k):
        check = self.columns >= k
        if check == False:
            raise ValueError ('index out of range')
        else:
            iden = self.Identidad(self.rows, self.columns)
            iden = iden.DelCol(k)            
            return self*iden
               
    
    def SwapRows(self,j,k):
        check = self.rows >= j and self.columns >= k
        if check == False:
            raise ValueError('rows or columns out of range')
        else:
            iden =self.Identidad(self.rows,self.columns)
            
            pos_1 = iden.GetPos(j, j)
            pos_2 = iden.GetPos(k, k)
            pos_3 = iden.GetPos(j, k)
            pos_4 = iden.GetPos(k, j)
            
            iden.elems[pos_1] = 0
            iden.elems[pos_2] = 0
            iden.elems[pos_3] = 1
            iden.elems[pos_4] = 1
            
            return (iden*self)
    
    def SwapCols(self,j,k):
        check = self.rows >= j and self.columns >= k
        if check == False:
            raise ValueError('rows or columns out of range')
        
        else:
            iden =self.Identidad(self.rows,self.columns)
        
            pos_1 = iden.GetPos(j, j)
            pos_2 = iden.GetPos(k, k)
            pos_3 = iden.GetPos(j, k)
            pos_4 = iden.GetPos(k, j)
            
            iden.elems[pos_1] = 0
            iden.elems[pos_2] = 0
            iden.elems[pos_3] = 1
            iden.elems[pos_4] = 1
        
            return (self*iden)

    def ScaleRow(self,j,x):
        new_mat = []
        for i in range(len(self.elems)):
            x_coord = self.GetCoords(i)[0]
            if x_coord != j:
                new_mat.append(self.elems[i])
            else:
                new_mat.append(self.elems[i]*x)
        return Myarray(new_mat, self.rows, self.columns, self.orden)
        
    def ScaleCol(self,k,y):
        new_mat = []
        for i in range(len(self.elems)):
            y_coord = self.GetCoords(i)[1]
            if y_coord != k:
                new_mat.append(self.elems[i])
            else:
                new_mat.append(self.elems[i]*y)
        return Myarray(new_mat, self.rows, self.columns)
                
    def Traspose(self):
        salida = Myarray([0]*self.rows*self.columns, self.columns,self.rows,self.orden)
        for i in range(len(self.elems)):
            coords = salida.GetCoords(i)
            salida.elems[i] = self.GetElem(coords[1],coords[0])
        return salida 
       
    def FlipRows(self):
        new_elems = []
        for i in range(self.rows,0,-1):
            row = self.GetRow(i)
            new_elems.extend(row)
        return Myarray(new_elems, self.rows, self.columns)
    
    def FlipCols(self):
        new_elems = []
        for i in range(self.columns,0,-1):
            col = self.GetCol(i)
            new_elems.extend(col)
        return Myarray(new_elems, self.rows, self.columns)
            
    
    def Det(self):
        if not self.rows == self.columns:
            raise ValueError('the number of rows must be equal to the number of columns')
        else:
            if self.rows == 1:
                return self.elems[0]            
            sumatoria = 0 
            for col in range(1, self.columns + 1):
                f_arbitraria = 1
                submatriz = self.GetSubmatrix([f_arbitraria],[col])
                sumatoria += self.GetElem(f_arbitraria,col)*(-1)**((f_arbitraria)+col)* submatriz.Det()                
            return sumatoria
    
    def Inversa(self):
        if not self.rows == self.columns:
            raise ValueError('the number of rows must be equal to the number of columns')
        else:
            if self.Det() == 0:
                raise ValueError ('this array cannot be inversed, det = 0')
            else:
                lista_identidades = []
                
                for col in range(1, self.columns + 1):
                    identidad = self.Identidad(self.rows, self.columns)
                    columna = self.GetCol(col)
                    pivote = self.GetElem(col, col)                    
                    
                    if pivote == 0:
                        index = col
                        for elem in columna[col::]:
                            if elem != 0:
                                pivote = elem  
                                index += 1
                                break
                            index += 1
                            
                        self = self.SwapRows(col, index)
                        permutador = self.Identidad(self.columns, self.rows).SwapRows(col, index)
                        lista_identidades.append(permutador)                    
                        columna = self.GetCol(col)
                        pivote = self.GetElem(col, col) 
                        
                    for i in range(len(columna)): 
                        columna[i] = -columna[i]
                    columna[col-1] = 1
                    for i in range(len(columna)):
                        columna[i] = columna[i]/pivote
                    
                    contador = 0
                    for i in range(len(identidad.elems)):    
                        if identidad.GetCoords(i)[1] == col:                    
                            identidad.elems[i] = columna[contador]
                            contador += 1
                    
                    lista_identidades.append(identidad)
                    self = identidad*self  
                    
                inversa = 1
                for iden in lista_identidades[::-1]:
                    inversa = inversa*iden
                
                return inversa
                    
    def Inversa_SinDet(self):
        if not self.rows == self.columns:
            raise ValueError('the number of rows must be equal to the number of columns')
        else:
            
            lista_identidades = []
            
            for col in range(1, self.columns + 1):
                identidad = self.Identidad(self.rows, self.columns)
                columna = self.GetCol(col)
                pivote = self.GetElem(col, col)                    
                
                if pivote == 0:
                    index = col
                    for elem in columna[col::]:
                        if elem != 0:
                            pivote = elem  
                            index += 1
                            break
                        index += 1
                        
                    self = self.SwapRows(col, index)
                    permutador = self.Identidad(self.columns, self.rows).SwapRows(col, index)
                    lista_identidades.append(permutador)                    
                    columna = self.GetCol(col)
                    pivote = self.GetElem(col, col) 
                    
                for i in range(len(columna)): 
                    columna[i] = -columna[i]
                columna[col-1] = 1
                for i in range(len(columna)):
                    columna[i] = columna[i]/pivote
                
                contador = 0
                for i in range(len(identidad.elems)):    
                    if identidad.GetCoords(i)[1] == col:                    
                        identidad.elems[i] = columna[contador]
                        contador += 1
                
                lista_identidades.append(identidad)
                self = identidad*self  
                
            inversa = 1
            for iden in lista_identidades[::-1]:
                inversa = inversa*iden
            
            return inversa
        
    def factor_LU(self):
        if self.rows != self.columns:
            raise ValueError("Matrix must be square for LU decomposition.")
          
        L = self.Identidad(self.rows, self.columns)
        U = Myarray(list(self.elems), self.rows, self.columns)      
        
        for i in range(self.rows):            
            for j in range(i+1, self.rows):
                factor = U.GetElem(j+1, i+1) / U.GetElem(i+1, i+1)
                L.elems[L.GetPos(j+1, i+1)] = factor  
                for k in range(1, self.columns+1):      
                    U.elems[U.GetPos(j+1, k)] -= factor * U.GetElem(i+1, k)
        
        return L, U


class Vector(Myarray):
    def __init__(self, elems):
        super().__init__(elems, len(elems), 1)
        
    def inner(self, other):
        return (self.Traspose() * other.Traspose().Traspose()).elems[0]
    
    def outer(self, other):
        return self * other.Traspose()
        
    def norm_2(self):
        return self.inner(self)
            
    def norma(self):
        return ((self.norm_2()**0.5))
    
    def versor(self):
        if self.norma() == 0:
            raise ValueError('la norma es 0')
        versor = self*(1/self.norma())
        return versor
        
    def proyeccion(self, other):
        versor = other.versor()
        b = self.inner(versor)
        resultado = b*versor
        return b, versor, resultado
    
    def orthogonal(self, other):
        proyecciones = self.proyeccion(other)
        proyeccion = proyecciones[2]
        return self - proyeccion  
    
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            resultado = super().__mul__(other)
            return Vector(resultado.elems)
        else:
            return super().__mul__(other)
    
    def __add__(self, other):
        return Vector((super().__add__(other).elems))

    def __sub__(self, other):
        return Vector((super().__sub__(other).elems))
    
    
class opt2d():
    def __init__(self, f, hx=0.001, hy=0.001):
        self.f = f
        self.hx = hx
        self.hy = hy
    
    def fx(self, x):
        dfx = (self.f(x[0] + self.hx, x[1]) - self.f(x[0] - self.hx, x[1])) / (2 * self.hx)
        return dfx
    
    def fy(self, x):
        dfy = (self.f(x[0], x[1] + self.hy) - self.f(x[0], x[1] - self.hy)) / (2 * self.hy)
        return dfy
    
    def fxx(self, x):
        dfxx = (self.fx((x[0] + self.hx, x[1])) - self.fx((x[0] - self.hx, x[1]))) / (2 * self.hx)
        return dfxx
    
    def fyy(self, x):
        dfyy = (self.fy((x[0], x[1] + self.hy)) - self.fy((x[0], x[1] - self.hy))) / (2 * self.hy)
        return dfyy
    
    def fxy(self, x):
        dfxy = (self.f(x[0] + self.hx, x[1] + self.hy) - self.f(x[0] + self.hx, x[1] - self.hy) - self.f(x[0] - self.hx, x[1] + self.hy) + self.f(x[0] - self.hx, x[1] - self.hy)) / (4 * self.hx * self.hy)
        return dfxy
        
    def gradf(self, x):
        dfx = self.fx(x)
        dfy = self.fy(x)
        return (dfx, dfy)

    def fv(self,v,x):
        vector_v = Vector(v)
        vector_gradiente = Vector(self.gradf(x))
        derivada_direccional = vector_gradiente.inner(vector_v.versor())
        return derivada_direccional
      
    def campo_gradiente(self, x_range, y_range, nx, ny, scale = 0.2, **kwargs):
        x_values = [x_range[0] + i *(x_range[1] - x_range[0]) / (nx - 1) for i in range(nx)]
        y_values = [y_range[0] + i *(y_range[1] - y_range[0]) / (ny - 1) for i in range(ny)]
        for x in x_values:
            for y in y_values:
                grad = self.gradf((x,y))
                plt.arrow(x, y, scale*grad[0], scale*grad[1], **kwargs)
        
        plt.xlim(x_range)
        plt.ylim(y_range)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Campo de Gradiente')
        plt.grid()
        

    
    def contour1(self, x0, y0, x_range):
        self.k = self.f(x0, y0)
        points = [(x0, y0)]
        xi_values = np.linspace(x_range[0], x_range[1], 500)
        self.x = x0
        self.y = y0

        def func(y):
            return self.f(self.x, y) - self.k

        def func_prime(y):
            return self.fy((self.x , y))

        for x in xi_values:
            self.x = x
            newton_solver = NewtonMethod()
            y_sol = newton_solver.find_root(func, func_prime, x0=self.y)
            if y_sol is None:
                continue
            self.y = y_sol
            points.append((self.x, self.y))
            if abs(func(-y_sol)) < 0.0001:
                points.append((self.x, -self.y))

        x_values, y_values = zip(*points)
        plt.scatter(x_values, y_values, label=f'k = {self.k}')
        plt.scatter(x0, y0, color='red')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Curva de nivel')
        plt.legend()
        plt.grid()
        plt.show()
        
    
    def contour2(self, x0 , y0 , x_limits, y_limits, alpha = 0.01, n_max = 5000, tol = 0.001):
        x1 = x0
        y1 = y0
        history = [(x1,y1)]
        
        for n in range(n_max):
            dfx, dfy = self.gradf((x1, y1))
            v = Vector((-dfy, dfx))
            v_versor = v.versor()
            x1 = x1 + alpha * v_versor.elems[0]
            y1 = y1 + alpha * v_versor.elems[1]
            history.append((x1,y1))
            if (x1 - x0)**2 + (y1 - y0)**2 < tol**2:
                break
        
        x_values, y_values = zip(*history)
        plt.plot(x_values, y_values)
        
        x1 = x0
        y1 = y0
        history = [(x1,y1)]   
        for n in range(n_max):
            dfx, dfy = self.gradf((x1, y1))
            v = Vector((dfy, -dfx))
            v_versor = v.versor()
            x1 = x1 + alpha * v_versor.elems[0]
            y1 = y1 + alpha * v_versor.elems[1]
            history.append((x1,y1))
            if (x1 - x0)**2 + (y1 - y0)**2 < tol**2:
                break
        x_values, y_values = zip(*history)
        plt.plot(x_values, y_values, color='red')
        
        plt.scatter([x0], [y0], color='red')  
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Curva de Nivel')
        plt.xlim(x_limits)
        plt.ylim(y_limits)
        plt.grid(True)
    
    def gdescent(self, x0, y0, delta = 0.1, n_max = 1000, tol = 0.00001):
        x, y = x0, y0
        history = [(x, y)]
        for n in range(n_max):
            dfx,dfy = self.gradf((x, y))
            x1 = x - delta * dfx
            y1 = y - delta * dfy 
            distancia = (((x1-x)**2 + (y1-y)**2)**0.5)
            x,y = x1,y1
            history.append((x, y))
            if distancia < tol:
                break                  
        return (x,y), history
    
    def Hessiano(self, punto, decimals=5):
        fxx = self.fxx(punto)
        fxy = self.fxy(punto)
        fyx = fxy
        fyy = self.fyy(punto)
        
        matriz_hessiana = Myarray([fxx,fxy,fyx,fyy], 2, 2)
        matriz_hessiana.elems = [round(elem, decimals) for elem in matriz_hessiana.elems]
        
        return matriz_hessiana
    
    
    def polytopes(self):
        
        def check_points(puntos, funcion_evaluadas, tol = 0.0000001):
            x1 = Myarray(puntos[0], 1, 2)
            x2 = Myarray(puntos[1], 1, 2)
            x3 = Myarray(puntos[2], 1, 2)
            
            v1 = x2 - x1
            v2 = x3 - x1
            elems = v1.elems + v2.elems
            mat = Myarray(elems, 2, 2)
            is_valid = True
            
            if mat.Det() == 0:
                is_valid = False
            
            distancia_x1_x2 = ( (x1.elems[0] - x2.elems[0])**2 + (x1.elems[1] - x2.elems[1])**2 ) **(1/2) 
            distancia_x2_x3 = ( (x2.elems[0] - x3.elems[0])**2 + (x2.elems[1] - x3.elems[1])**2 ) **(1/2) 
            distancia_x1_x3 = ( (x1.elems[0] - x3.elems[0])**2 + (x1.elems[1] - x3.elems[1])**2 ) **(1/2) 
            if max(distancia_x1_x2, distancia_x2_x3, distancia_x1_x3) < tol:
                is_valid = False
            
            fa, fb, fc = funcion_evaluadas
            fa_fb = abs(fa - fb)
            fb_fc = abs(fb - fc)
            fa_fc = abs(fa - fc)
            if max(fa_fb, fb_fc, fa_fc) < tol:
                is_valid = False
            
            return is_valid 
        
        
        def reflexion(puntos):
            v1 = Vector(puntos[-1]) - Vector(puntos[-2])
            versor_v1 = v1.versor()
            
            aux = Vector(puntos[-1]) - Vector(puntos[0])
            a = aux.inner(versor_v1)
            b = aux - a*versor_v1
            
            reflexion = Myarray(puntos[0], 1, 2) + 2*b
            
            return tuple(reflexion.elems)
        
        def shrink(puntos):
            x1_prime = (0.5*(puntos[0][0] + puntos[-1][0]), 0.5*(puntos[0][1] + puntos[-1][1]))
            x2_prime = (0.5*(puntos[1][0] + puntos[-1][0]), 0.5*(puntos[1][1] + puntos[-1][1]))
            new_points = [x1_prime,x2_prime,puntos[-1]]
            sorted(new_points)
            
            return new_points
        
        #puntos = [(1,0), (2,1), (3,4)]
        puntos = [(random.uniform(-5, 5), random.uniform(-5, 5)) for _ in range(3)]
        puntos_ordenados = puntos
        history = [(puntos)]
        
        for i in range(1000):
            
            puntos_ordenados = sorted(puntos_ordenados, key=lambda p: self.f(*p), reverse=True)
            
            fa = self.f(puntos_ordenados[0][0], puntos_ordenados[0][1])
            fb = self.f(puntos_ordenados[1][0], puntos_ordenados[1][1])
            fc = self.f(puntos_ordenados[2][0], puntos_ordenados[2][1])
            check = check_points(puntos_ordenados, (fa,fb,fc))
            
            if not check:
                break 
            
            reflexion_x = reflexion(puntos_ordenados)
            
            if self.f(reflexion_x[0],reflexion_x[1]) < self.f(puntos_ordenados[0][0], puntos_ordenados[0][1]):
                puntos_ordenados[0] = reflexion_x
                
            else: 
                puntos_switch = [puntos_ordenados[1], puntos_ordenados[0], puntos_ordenados[-1]]
                reflexion_x2 = reflexion(puntos_switch)
                if self.f(reflexion_x2[0],reflexion_x2[1]) < self.f(puntos_ordenados[1][0], puntos_ordenados[1][1]):
                    puntos_ordenados[1] = reflexion_x2
                
                else: 
                    puntos_ordenados = shrink(puntos_ordenados)
            
            history.append(list(puntos_ordenados))
            
            optimo = (puntos_ordenados[-1][0], puntos_ordenados[-1][1])
            
        plt.figure()
        x_range = (optimo[0] - 5, optimo[0] + 5)
        y_range = (optimo[1] - 5, optimo[1] + 5)
        
        
        self.campo_gradiente(x_range, y_range, 20, 20)
        
        self.contour2(optimo[0], optimo[1], x_range, y_range) 
        
        for puntos in history:
            x_values, y_values = zip(*puntos)
            plt.plot(x_values, y_values, 'o-')
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Evolución de los puntos en el método de polytopes')
        plt.grid()
        plt.show()
            
        return (round(optimo[0], 2), round(optimo[1], 2)), history

        
class opt2d_con_restriccion(opt2d):
    def __init__(self, f, restriccion, hx = 0.001, hy = 0.001):
        super().__init__(f, hx, hy)
        self.restriccion = restriccion
    
    
    def minimize(self, punto, method = 'polytopes'):
        
        def busqueda_sobre_restriccion(punto):
            func_restriccion_cuadrada = lambda x,y: self.restriccion(x,y)**2 
            func_restriccion_cuadrada_opt = opt2d(func_restriccion_cuadrada)
            if method == 'polytopes':
                optimo, history = func_restriccion_cuadrada_opt.polytopes()
            else:
                optimo, history = func_restriccion_cuadrada_opt.gdescent(punto[0], punto[1])
            return optimo, history
        
        
        def navegar_restriccion(alpha = 0.01):
            xn, history = busqueda_sobre_restriccion(punto)
            restriccion = opt2d(self.restriccion)
            history_res = [xn]
            for _ in range(1000):
                grad_fx = Vector(self.gradf(xn))
                grad_gx = Vector(restriccion.gradf(xn)).versor()
                a = (grad_gx.inner(-1*grad_fx)) * grad_gx
                vec = -1*grad_fx - a
                
                xn1 = (xn[0] + alpha * vec.elems[0], xn[1] + alpha * vec.elems[1])
                if  self.f(xn1[0], xn1[1]) > self.f(xn[0], xn[1]):
                    history_res.append(xn)
                    break 
                history_res.append(xn1)
                xn = xn1
            return xn, history, history_res
    
        return navegar_restriccion()
    
class opt2d_restriccion_desigualdad(opt2d):
    def __init__(self, f, desigualdad, hx = 0.001, hy = 0.001):
        super().__init__(f, hx, hy)
        self.desigualdad = desigualdad
 
    def minimize(self, punto, method = 'polytopes'):   
         
        def busqueda_zona_admisible(punto):
            func_desigualdad = lambda x,y: min(self.desigualdad(x,y),0)**2 
            func_restriccion_cuadrada_opt = opt2d(func_desigualdad)
            if method == 'polytopes':
                optim, history = func_restriccion_cuadrada_opt.polytopes()
            else:
                optim, history = func_restriccion_cuadrada_opt.gdescent(punto[0], punto[1])
            return optim, history
            
            
        def navegar_restriccion(alpha = 0.01):
            xn, history = busqueda_zona_admisible(punto)
            restriccion = opt2d(self.desigualdad)
            history_res = [xn]
            for _ in range(1000):
                grad_fx = Vector(self.gradf(xn))
                grad_gx = Vector(restriccion.gradf(xn)).versor()
                a = (grad_gx.inner(-1*grad_fx)) * grad_gx
                vec = -1*grad_fx - a
                
                xn1 = (xn[0] + alpha * vec.elems[0], xn[1] + alpha * vec.elems[1])
                if  self.f(xn1[0], xn1[1]) > self.f(xn[0], xn[1]):
                    history_res.append(xn)
                    break 
                history_res.append(xn1)
                xn = xn1
            return xn, history, history_res
        
        return navegar_restriccion(self)
    
    
    
    
    
func = lambda x,y: y 
    
res = lambda x,y : (x-1)**2 + (y-2)**2 - 1 
#res = lambda x,y: y**3 - x**2 

opt = opt2d_restriccion_desigualdad(func,res)
opt.minimize((4,2), 'gdesdcent')








