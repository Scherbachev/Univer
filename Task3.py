# -*- coding: utf-8 -*-

import numpy
import tools

class GaussianPlaneWave:
    ''' Класс с уравнением плоской волны для гауссова сигнала в дискретном виде
    d - определяет задержку сигнала.
    w - определяет ширину сигнала.
    Sc - число Куранта.
    eps - относительная диэлектрическая проницаемость среды, в которой расположен источник.
    mu - относительная магнитная проницаемость среды, в которой расположен источник.
    '''
    def __init__(self, d, w, Sc=1.0, eps=1.0, mu=1.0):
        self.d = d
        self.w = w
        self.Sc = Sc
        self.eps = eps
        self.mu = mu
        
    def getE(self, m, q):
        '''
        Расчет поля E в дискретной точке пространства m
        в дискретный момент времени q
        '''
        return numpy.exp(-(((q - m * numpy.sqrt(self.eps * self.mu) / self.Sc) - self.d) / self.w) ** 2)

if __name__ == '__main__':
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * numpy.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 800

    # Размер области моделирования в м
    X = 5.5

    # Размер ячейки по пространству
    dx = 0.01

    # Размер области моделирования в отсчетах
    maxSize = int(X/dx)

    # Шаг по времени
    dt = Sc * dx / 300000000
    
    # Датчики для регистрации поля
    probesPos = [int(maxSize*0.75)]
    probes = [tools.Probe(pos, maxTime) for pos in probesPos]

    # Положение источника
    sourcePos = int(maxSize/2)

    # Параметры среды
    eps = numpy.ones(maxSize)
    eps = 3.5 * eps
    mu = numpy.ones(maxSize-1)

    # Начало идеально поглощающего слоя
    layer_loss_x = (maxSize - 100)

    # Потери в среде
    loss = numpy.zeros(maxSize)
    loss[layer_loss_x:] = 0.02

    # Коэффициенты для расчета поля E
    ceze = (1 - loss) / (1 + loss)
    cezh =  Sc * W0 / (eps * (1 + loss))

    # Коэффициенты для расчета поля H
    chyh = (1 - loss) / (1 + loss)
    chye = Sc / (W0 * (1 + loss))

    # Усреднение коэффициентов на границе поглощающего слоя
    ceze[layer_loss_x] = (ceze[layer_loss_x - 1] + ceze[layer_loss_x + 1]) / 2
    cezh[layer_loss_x] = (cezh[layer_loss_x - 1] + cezh[layer_loss_x + 1]) / 2
    
    Ez = numpy.zeros(maxSize)
    Hy = numpy.zeros(maxSize-1)

    source = GaussianPlaneWave(45, 20, Sc, eps[sourcePos], mu[sourcePos])

    for probe in probes:
        probe.addData(Ez, Hy)

    # Создание экземпляра класса для отображения
    # распределения поля в пространстве
    display = tools.AnimateFieldDisplay(maxSize, -1.1, 1.1, 'Ez, В/м', dx)
    display.activate()
    display.drawSources([sourcePos])
    display.drawProbes(probesPos)

    for q in range(1,maxTime):
        # Расчет компоненты поля H
        Hy =  Hy * chyh[:-1] + chye[:-1] * (Ez[1:] - Ez[:-1]) 

        Hy[sourcePos-1] -= Sc / (W0 * mu[sourcePos-1]) * source.getE(0, q)

        # Источник возбуждения
        Ez[sourcePos] += Sc / numpy.sqrt(eps[sourcePos]*mu[sourcePos]) * source.getE((1/eps[sourcePos])-1, q+(1/eps[sourcePos])**2)
        
        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:-1] = ceze[1: -1] * Ez[1:-1] + cezh[1: -1] * (Hy[1:] - Hy_shift)


        # Регистрация поля в датчиках
        for probe in probes:
            probe.addData(Ez, Hy)

        if q % 2 == 0:
            display.updateData(Ez, q)

    display.stop()

    # Отображение сигнала, сохраненного в датчиках
    tools.showProbeSignals(probes, -1.1, 1.1, dt)

    # Отображение спектра сигнала
    spectr = tools.Spectrum(probe.E, dt, 3e9)
    spectr.fourierTransform()
    spectr.ShowNorm()
