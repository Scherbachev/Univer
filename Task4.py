# -*- coding: utf-8 -*-
import pylab
import numpy
import tools
import math

class GaussianPlaneWaveM:
    ''' Класс с уравнением плоской волны для гауссова сигнала в дискретном виде
    d - определяет задержку сигнала.
    w - определяет ширину сигнала.
    Sc - число Куранта.
    eps - относительная диэлектрическая проницаемость среды, в которой расположен источник.
    mu - относительная магнитная проницаемость среды, в которой расположен источник.
    '''
    def __init__(self, d, w, Tm,  Sc=1.0, eps=1.0, mu=1.0):
        self.d = d
        self.w = w
        self.Tm = Tm
        self.Sc = Sc
        self.eps = eps
        self.mu = mu
        
    def getE(self, m, q):
        '''
        Расчет поля E в дискретной точке пространства m
        в дискретный момент времени q
        '''
        return numpy.sin(2*math.pi*q/Tm) * numpy.exp(-(((q - m * numpy.sqrt(self.eps * self.mu) / self.Sc) - self.d) / self.w) ** 2)

if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 6000

    # Размер области моделирования в м
    X = 0.5

    # Размер ячейки по пространству
    dx = 0.00025

    # Размер области моделирования в отсчетах
    maxSize = int(X/dx)

    # Шаг по времени
    dt = Sc * dx / 300000000
    
    # Датчики для регистрации поля
    probesPos = [25,75]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Положение источника
    sourcePos = int(50)

    #Свойства слоев диэлектрика
    d1 = int(0.04/dx)
    d2 = int(0.08/dx)
    eps1 = 1.5
    eps2 = 5.9
    eps3 = 2.6

    #Координаты начала слоев
    lay_st = int(maxSize/2)
    lay_1 = lay_st + d1
    lay_2 = lay_1 + d2
    
    # Параметры среды
    eps = numpy.ones(maxSize)
    eps[lay_st:lay_1] = eps1
    eps[lay_1:lay_2] = eps2
    eps[lay_2:] = eps3
    mu = numpy.ones(maxSize-1)

    # Коэффициенты для расчета E и H полей
    K_E = Sc * W0 / eps
    K_H = Sc / (W0 * mu)

    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize-1)

    # Создание источника поля
    T = 60 # Половина длительности импульса (В отчетах)
    Tm = 120 # Период несущей (В отчетах)
    source = GaussianPlaneWaveM(2*T + 20, T, Tm, Sc, eps[sourcePos], mu[sourcePos])

    for probe in probes:
        probe.addData(Ez, Hy)

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(maxSize, -1.1, 1.1, 'Ez, В/м', dx)
    display.activate()
    display.drawSources([sourcePos])
    display.drawProbes(probesPos)
    display.drawBoundary(lay_st)
    display.drawBoundary(lay_1)
    display.drawBoundary(lay_2)

    # Расчет коэффициентов для граничных условий
    K_L = (((Sc/numpy.sqrt(eps[0]*mu[0]))-1)/((Sc/numpy.sqrt(eps[0]*mu[0]))+1))
    K_R = (((Sc/numpy.sqrt(eps[-1]*mu[-1]))-1)/((Sc/numpy.sqrt(eps[-1]*mu[-1]))+1))

    #
    Ez_oldL = Ez[1]
    Ez_oldR = Ez[-1]

    for q in range(1,maxTime):
        # Расчет компоненты поля H
        Hy =  Hy  + K_H * (Ez[1:] - Ez[:-1]) 

        Hy[sourcePos-1] -= Sc / (W0 * mu[sourcePos-1]) * source.getE(0, q)

        # Источник возбуждения
        Ez[sourcePos] += Sc / numpy.sqrt(eps[sourcePos]*mu[sourcePos]) * source.getE(-0.5, q+0.5)
        
        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = Ez[1:-1] + K_E[1:-1] * (Hy[1:] - Hy_shift)
        Ez[0] = Ez_oldL + K_L * (Ez[1] - Ez[0])
        Ez_oldL = Ez[1]
        Ez[-1] = Ez_oldR + K_R * (Ez[-2] - Ez[-1])
        Ez_oldR = Ez[-2]

        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 20 == 0:
            display.updateData(Ez, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    tools.showProbeSignals(probes, -1.1, 1.1, dt)

    # Исключение отраженных волн из датчика падающей волны с помощью временного окна
    # Временное окно равно длительности импульса + времени, необходимому волне для прохождения
    # расстояния от левого края до диэлектрика (в вакууме)
    probes[1].E[int(T + lay_st*dx/(300000000*dt)):] = 0

    # Формирование спектра падающей и отраженной волн
    spectr1 = tools.Spectrum(probes[0].E, dt, 30e9)
    spectr1.fourierTransform()
    spectr2 = tools.Spectrum(probes[1].E, dt, 30e9)
    spectr2.fourierTransform()

    # Вывод спектров падающей и отраженной волн, нормированных к максимуму спектра падающей волны
    max_Norm = max(spectr2.PF)
    fig, ax = pylab.subplots()
    ax.plot(spectr1.f, spectr1.PF/max_Norm)
    ax.plot(spectr1.f, spectr2.PF/max_Norm)
    ax.set_xlim(0, spectr1.xMax)
    ax.set_xlabel('f, Гц')
    ax.set_ylabel(r'|$\frac{P}{P_{max}}$|')
    ax.grid()
    pylab.show() 

    # Вывод модуля коэффициента отражения в полосе частот Fmin, Fmax
    Fmin = 5e9
    Fmax = 15e9
    fig, ax = pylab.subplots()
    ax.plot(spectr1.f, spectr1.PF/spectr2.PF)
    ax.set_xlim(Fmin, Fmax)
    ax.set_ylim(0, 1)
    ax.set_xlabel('f, Гц')
    ax.set_ylabel('|Г|')
    ax.grid()
    pylab.show() 
