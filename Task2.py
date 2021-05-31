import math
from scipy.special import spherical_jn, spherical_yn
import numpy as np
import requests as rqst
import matplotlib.pyplot as plt
import os
import re

# Рассчитать ЭПР
# Построить график
# Сохранить результаты в файл
# Radar cross section

#Функция, которая из строк вида Параметр="MeP" и Параметр="X" возвращает число
def desif(S):
    S = S[S.find('=')+2:-1] #Вырезаем название параметра и кавычки
    if S.find('e') == -1:   #Если в кавычках число не представлено в инженерном формате
        return float(S)     #То мы возвращаем число
    else:
        S = S.split('e')    #Иначе разделяем его на мантису и порядок
        return float(S[0]) * (10 ** float(S[1]))    #Переводим в десятичный формат и возвращаем его

#Финкция, которая из файла по ссылке получает исходные данные
#Считаем, что файл доступен и в нем гарантированно есть одна строка с нужным вариантом
#(Не осуществляем проверку на такие сценарии)
def request_www(n ,link):
    file = link.split('/')[-1]  #Имя файла вытаскиваем из последней части ссылки, после последнего слеша
    r = rqst.get(link, allow_redirects = True)  #Скачиваем файл
    open(file, 'wb').write(r.content)               #Записываем файл
    with open(file) as f:
        result = []                                 #Создаем пустой массив под результат
        for l in f:                                 #Пробегаемся по всем строкам
            if l.find('number="'+str(n)+'"') != -1: #Если на строке указан номер варианта, то
                lt = l[2:-3].split(' ')             #Разбиваем её на отдельные параметры(За вычетом нач. пробелов и символов в конце)
                result.append(desif(lt[2]))         #Записываем их в массив по фиксированным индесам
                result.append(desif(lt[3]))
                result.append(desif(lt[4]))
                return result

#Функция для вычисления ЭПР
def cntRCS(lam, r):
    summ = 0
    kr = 2 * math.pi * r / lam
    # Задаем значения функций Бесселя для n = 0 для первой итерации
    J_prev = spherical_jn (0, kr)
    Y_prev = spherical_yn (0, kr)
    H_prev = J_prev + 1j * Y_prev
    for n in range(1, 50):
        # Вычисляем значения функций Бесселя для текущей n
        J_now = spherical_jn (n, kr)
        #J_prev = spherical_jn (n - 1, kr)
        Y_now = spherical_yn (n, kr)
        #Y_prev = spherical_yn (n - 1, kr)
        H_now = J_now + 1j * Y_now
        #H_prev = J_prev + 1j * Y_prev
        # Считаем коэффициенты a и b
        a = J_now / H_now
        b = (kr * J_prev - n * J_now) / (kr * H_prev - n * H_now)
        summ += ((-1) ** n) * (n + 0.5) * (b - a)
        # Переносим значения функций Бесселя на следующий шаг
        J_prev = J_now
        Y_prev = Y_now
        H_prev = H_now
    return lam * lam * np.abs(summ) * np.abs(summ) / math.pi

#Функция, которая строит график
def graph(lambd, rcs):
    plt.plot(lambd, rcs)
    plt.xlabel('wavelength, m')
    plt.ylabel('RCS, m^2')
    plt.grid()
    plt.show()

#Начальные параметры

n = 11 #Номер варианта
N = 100  #Число точек на графике

#Основная часть программы
url = 'https://jenyay.net/uploads/Student/Modelling/task_02.xml'

initpar = request_www(n, url)

lambd_min = 300000000 / initpar[2]  #Минимальная длина волны определяется максимальной частотой
lambd_max = 300000000 / initpar[1]  #Максимальная длина волны определяется минимальной частотой
lambd = np.arange(start=lambd_min, stop=lambd_max, step=(lambd_max-lambd_min)/N)
RCS = cntRCS(lambd, initpar[0])

path = os.path.dirname(__file__)
path = os.path.join("results.txt")
with open(path, 'w') as f:
    for i in range(len(RCS)):
        f.write('%9.7f    %9.7f\n' % (lambd[i], RCS[i]))

graph(lambd, RCS)
