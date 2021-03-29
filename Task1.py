import os
import matplotlib.pyplot as plt
import math

# Начало и конец интервала, шаг дискретизации
Xs = -15
Xf = 5
dX = 0.01

# Список точек для определения значения функции и пустой список значений
X = [i*dX+Xs for i in range(int((Xf-Xs)/dX))]
Y = []

# Заполнение списка значений по списку точек
for i in range(len(X)):
	Temp = X[i]
	Y.append(100*((abs(1-0.01*Temp*Temp))**0.5)+0.01*abs(Temp+10))

# Определение положения файла
dirname = os.path.dirname(__file__)
# Проверка существования папки и ей создание при отсутствии
filename = os.path.join(dirname, 'results')
if os.path.exists(filename)==False:
    os.mkdir(filename)
# Создание файла, начало записи
filename = os.path.join(filename, 'results.xml')
f = open(filename,'w')
# Шапка файла
f.write('<?xml version="1.1" encoding="UTF-8" ?>\n')
# Данные
f.write('<data>\n')
f.write('\t<xdata>\n')
for i in range(len(X)):
    f.write('\t\t<x>' + str(X[i]) + '</x>\n')
f.write('\t</xdata>\n')
f.write('\t<ydata>\n')
for i in range(len(Y)):
    f.write('\t\t<y>' + str(Y[i]) + '</y>\n')
f.write('\t</ydata>\n')
f.write('</data>')
# Закрытие файла
f.close()
#Построение графика
plt.plot(X,Y)
plt.show()

