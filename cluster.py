# класс клетки
import copy
import itertools
import numpy as np
import matplotlib.pyplot as plt

class Cell:
    def __init__(self, superCell, localPosition):
        """ при инициализации определяются только объект надклетки и локальная позиция в ней """
        self.superCell = superCell
        self.depth = self.superCell.depth + 1
        self.localPosition = localPosition

    def specify(self, CONSTANTS):
        """ задаются определенные свойства на основе параметров, передаваемых в CONSTANT:
        lenght - длина ребра клетки, measure - мера клетки, и globalPosition """
        self.lenght = self.superCell.lenght / CONSTANTS['DENOMINATOR']
        self.measure = self.lenght ** CONSTANTS['DIMENSION']
        globalPosition = []
        for i in range(CONSTANTS['DIMENSION']):
            globalPosition.append(self.superCell.globalPosition[i] * CONSTANTS['DENOMINATOR'] + self.localPosition[i])
        self.globalPosition = tuple(globalPosition)

    def disperse(self, X, CONSTANTS):
        """ клетка разбивается на подклетки, при этом создаются словарь, содержащий подклетки
        и словарь, содержащий плотности клеток (по ключам, соответствующим локальной позиции"""
        self.subCells = {}
        self.subCellsDensity = {}
        s = list(range(CONSTANTS['DENOMINATOR']))

        subPositionList = list(itertools.product(s, repeat=CONSTANTS['DIMENSION']))
        sampleList = list(self.sampleList)
        for i in range(len(subPositionList)):
            self.subCells[subPositionList[i]] = Cell(self, subPositionList[i])
            sampleList = self.subCells[subPositionList[i]].filling(sampleList, X, CONSTANTS)
            self.subCellsDensity[subPositionList[i]] = self.subCells[subPositionList[i]].density
        return self.subCells

    def getBorders(self, CONSTANTS):
        """ возвращает матрицу границ клетки как массив numpy """
        borders = np.zeros((CONSTANTS['DIMENSION'], 2))
        for i in range(CONSTANTS['DIMENSION']):
            borders[i][0] = self.lenght * self.globalPosition[i]
            borders[i][1] = borders[i][0] + self.lenght
        return borders

    def filling(self, sampleListFree, X, CONSTANTS):
        ''' sampleListFree - список индексов доступных экземпляров
        X - матрица всех экземпляров
        словарь CONSTANTS содержит (SIZE, DIMENSION, DENOMINATOR)
        создает атрибут со списком экземпляров, вошедших в клетку (sampleList)
        возвращает список идексов экземпляров, не вошедших в клетку'''
        
        self.specify(CONSTANTS)

        sampleListRemain = copy.copy(sampleListFree)
        borders = self.getBorders(CONSTANTS)
        self.sampleList = []    # список с номерами попавших в клетку экземпляров
        for i in sampleListFree:
            for j in range(CONSTANTS['DIMENSION']):
                indicator = X[i][j] >= borders[j][0] and X[i][j] <= borders[j][1]
                if not indicator: break
            if indicator:
                sampleListRemain.remove(i)
                self.sampleList.append(i)
        self.sampleList = tuple(self.sampleList)
        self.density = len(self.sampleList) / self.measure
        return sampleListRemain
        
    def separation(self, depth, CONSTANTS):
        """ метод применяется для неоднородных клеток, если идет дальнейшее разбиение
        сильно разреженные клетки удаляются, однородные клетки переносятся
        в словарь homoSubCells, остальное остается в словаре subCells, которые
        возвращают методом в кортеже """
        self.homoCells = {}
        deadCells = []
        HOMOGENEITY, DENSITY = CONSTANTS["HOMOGENEITY"], CONSTANTS["DENSITY"]
        for key in self.subCells:
            if self.subCells[key].density < (self.density / DENSITY):     
                deadCells.append(key)
        for key in deadCells:
            del self.subCells[key]
        for key in self.subCells:
            # print(key, self.subCells[key].getHomogeneity(depth, CONSTANTS))
            if self.subCells[key].getHomogeneity(depth, CONSTANTS) > HOMOGENEITY:
                self.homoCells[key] = self.subCells[key]
        for key in self.homoCells:
            del self.subCells[key]
        return (self.subCells, self.homoCells)
    
    # def selection(self, depth, CONSTANTS):
    #     """ метод применяется для однородных клеток, если идет дальнейшее разбиение
    #     сильно разреженные подклетки удаляются из словаря subCells 
    #     формируется словарь граничных подклеток """
    #     deadCells = []
    #     self.innerCells = {}
    #     self.borderCells = {}
    #     self.depth = depth

    #     outerSubCells = self.getOuterSubCells(CONSTANTS)
    #     for key in outerSubCells:
    #         if outerSubCells[key].density < (self.density / 5):
    #             deadCells.append(key)
    #         else: self.innerCells[key] = self.outerSubCells[key]

    #     neighborList = self.getNeighbors()
    #     for key in deadCells:
    #         for k in neighborList:
    #             if k in self.innerCells:
    #                 self.borderCells[k] = self.subCells[k]

    #     for key in deadCells:
    #         del self.subCells[key]
    #     return self.borderCells

    def getOuterSubCells(self, CONSTANTS):
        outerSubCells = dict()
        for key in self.subCells:
            if (0 in key or CONSTANTS['DENOMINATOR']-1 in key):
                outerSubCells[key] = self.subCells[key]
        return outerSubCells

    def getGlobalPosition(self, CONSTANTS):
        globalPosition = []
        DENOMINATOR = CONSTANTS["DENOMINATOR"]
        for i in range(CONSTANTS["DIMENSION"]):
            globalPosition.append(self.superCell.globalPosition[i] * DENOMINATOR + self.localPosition[i])
        return tuple(globalPosition)
    def getDensityStd(self):
        return np.std(np.array(list(self.subCellsDensity.values())))
    def getHomogeneity(self, depth, CONSTANTS):
        s = np.std(np.array(list(self.subCellsDensity.values())))
        if s > 0:
            return self.density * depth * CONSTANTS['DENOMINATOR'] / s
        else:
            return 0
    def getDensityMdn(self):
        return np.median(np.array(list(self.subCellsDensity.values())))
    def getMinMaxDnst(self):
        keys = tuple(self.subCellsDensity.keys())
        values = np.array(list(self.subCellsDensity.values()))
        min = np.min(values)
        max = np.max(values)
        argmax = np.argmax(values)
        argmin = np.argmin(values)
        return ((keys[argmin], keys[argmax]), (min, max))

    def getMinMaxDiff(self):
        values = np.array(list(self.subCellsDensity.values()))
        return np.max(values) - np.min(values)
    
    def expand(self, cellDict):
        neighborList = self.getNeighbors()
        mergeList = []
        for key in neighborList:
            if key in cellDict.keys():
                mergeList.append(key)
        return mergeList

    def getNeighbors(self):
        """ окрестность Неймана (позиции клеток, смежных по гиперграням) """
        dimension = len(self.globalPosition)
        neighborNeumanList = list()
        globalPosition = np.array(list(self.globalPosition))
        for i in range(dimension):
            globalPosition[i] -= 1
            neighborNeumanList.append(tuple(globalPosition))
            globalPosition[i] += 2
            neighborNeumanList.append(tuple(globalPosition))
            globalPosition[i] -= 1
        return neighborNeumanList

    def getExtraNeighbors(self):
        """ позиции клеток, смежных по гиперребрам """
        dimension = len(self.globalPosition)
        neighborExtraList = list()
        globalPosition = np.array(list(self.globalPosition))
        for i in range(dimension-1):
            for j in range(dimension):
                if j>i:
                    globalPosition[i] -= 1
                    globalPosition[j] -= 1
                    neighborExtraList.append(tuple(globalPosition))
                    globalPosition[i] += 2
                    neighborExtraList.append(tuple(globalPosition))
                    globalPosition[j] += 2
                    neighborExtraList.append(tuple(globalPosition))
                    globalPosition[i] -= 2
                    neighborExtraList.append(tuple(globalPosition))
                    globalPosition[i] += 1
                    globalPosition[j] -= 1
            return neighborExtraList

    # def getVertexNeighbors(self):
    #     """ окрестности вершин """
    #     dimension = len(self.globalPosition)
    #     shiftList = np.array(list(itertools.product([1,-1], repeat=dimension)))
    #     vertexNeighborList = list()
    #     globalPosition = np.array(list(self.globalPosition))
    #     for i in range(2 ** dimension):
    #             vertexNeighborList.append(tuple(globalPosition + shiftList[i]))
    #     return vertexNeighborList


class HeadCell(Cell):
    def __init__(self, CONSTANTS):
        self.localPosition = tuple([0] * CONSTANTS['DIMENSION'])
        self.globalPosition = self.localPosition
        self.depth = 0

    def specify(self, CONSTANTS):
        """ определяются конеркетные поля
        словарь CONSTANTS содержит (SIZE, DIMENSION, DENOMINATOR) """
        self.lenght = CONSTANTS['SIZE']
        self.measure = self.lenght ** CONSTANTS['DIMENSION']


class Cluster():
    def __init__(self, cellList, depth, CONSTANTS):
        """ при инициализации кластера первым аргументом передается список слившихся клеток и глубина """
        self.passiveCellDict = dict() # словарь пассивных клеток, ключ соответствует глубине
        self.activeCellDict = dict()            # словарь активных клеток
        self.passiveCellDict[depth] = dict()    # словарь пассивных клеток глубины depth
        self.activeCellDict[depth] = dict()     # словарь активных клеток глубины depth
        cellPstns = self.getGlobalPstns(cellList) # список глобальных позиций всех клеток
        # распределим все клетки по словарям для текущей глубины 
        for cell in cellList:
            for n in cell.getNeighbors():
                pstn = cell.getGlobalPosition(CONSTANTS)
            # если хотя бы одна из соседних позиций отсутствует в списке клеток кластера:
                if not (n in cellPstns):
                    self.activeCellDict[depth][pstn] = cell # ключ будет соответствовать глубине 
                    break
                else:   # иначе клетка будет пассивной
                    self.passiveCellDict[depth][pstn] = cell

    def getGlobalPstns(self, cellList):
        """ Принимет список клеток и возвращает список глобальных позцици"""
        cellPstns = list()
        for cell in cellList:
            cellPstns.append(cell.globalPosition)
        return cellPstns

    def generateFacePstns(self, axis, CONSTANTS):
        """ генерируются локальные позиции лицевых подклеток клетки вдоль 
        axis расположенных по оси axis, выдается два списка позиций для 
        каждой из противоположных сторон """
        DENOMINATOR = CONSTANTS["DENOMINATOR"]
        DIMENSION = CONSTANTS["DIMENSION"]
        s = list(range(CONSTANTS["DENOMINATOR"])) # список позиций всех подклеток
        subCellPstn = list(itertools.product(s, repeat=DIMENSION))
        faceCellPstnsSide_1 = list()    # позиции подклеток "левой" стороны
        faceCellPstnsSide_2 = list()    # позиции подклеток "правой" стороны
        for pstn in subCellPstn:        # отбираются подклетки сторон
            if pstn[axis] == 0:         # для "левой" стороны
                faceCellPstnsSide_1.append(pstn)
            elif pstn[axis] == DENOMINATOR - 1:    # для "правой" стороны
                faceCellPstnsSide_2.append(pstn)
        return (faceCellPstnsSide_1, faceCellPstnsSide_2)

    def collectOuterSubCells(self, depth, X, CONSTANTS):
        """ после того как кластер инициализирован, собираются все внешние подклетки
        активных клеток глубины depth; эти подклетки будут включены в словарь
        activeCellDict по ключу depth+1 """

        outerSubCells = dict()         # список, в котором соберутся активные подклетки
        # цикл перебирает все активные клетки уровня depth:
        for cell in self.activeCellDict[depth].values():
            neighborList = cell.getNeighbors()
            cell.disperse(X, CONSTANTS)         # дробление клетки
            
            # собираем список глобальных позиций клеток уровня depth:
            cellPstns = self.getGlobalPstns(list(self.passiveCellDict[depth].values()))
            cellPstns += self.getGlobalPstns(list(self.activeCellDict[depth].values()))

            outerSubCellPstns = list()  # в этом списке соберутся все локальные позиции 
            # внешних подклеток из данной клетки; для каждой оси axis собираем локальные 
            # позиции подклеток соответствующих сторон, если они внешние
            for axis in range(CONSTANTS['DIMENSION']):
                if not (neighborList[axis*2] in cellPstns):
                    outerSubCellPstns += self.generateFacePstns(axis, CONSTANTS)[0]
                if not (neighborList[axis*2+1] in cellPstns):
                    outerSubCellPstns += self.generateFacePstns(axis, CONSTANTS)[1]
            
            # имея в outerSubCellPstns локальные позиции внешних подклеток, 
            # собираем внешние подклетки в словарь по ключу глобальной позиции:
            for pstn in outerSubCellPstns:
                globalPosition = cell.subCells[pstn].getGlobalPosition(CONSTANTS)
                outerSubCells[globalPosition] = cell.subCells[pstn]
            
            # изменение словарей activeCellDict и passiveCellDict:
            self.activeCellDict[depth+1] = outerSubCells
            # добавляется пустой словарь по ключу depth+1, если такой ключ отсутствует
            self.passiveCellDict.setdefault(depth+1, dict())
        return outerSubCells
    
    def expand(self, remainHomoCells, otherCells, depth, CONSTANTS):
        """ метод расширения кластера на один шаг на глубине depth
        получает словарь оставшихся гомогенных клеток remainHomoCells, из которого 
        извлекутся все клетки, которые примкнут к данному кластеру, вторым
        аргументом принимается словарь негомогенных непустых клеток otherCells """
        
        # получим позиции активных клеток
        keys = list(self.activeCellDict[depth].keys())
        newCells = dict()       # в словарь newCells будут собираться клетки

        for key in keys:
            isFrozen = True     # индикатор замороженности клетки
            # проверим соседей, смежных по гиперребрам и если есть подходящие,
            # то они добавляются в словарь newCells:
            extraNgbrs = self.activeCellDict[depth][key].getExtraNeighbors()
            for pstn in extraNgbrs:
                if pstn in remainHomoCells:
                    adjacents = getAdjacentPair(key, pstn)
                    threshold = CONSTANTS["ADJACENCY"] * self.activeCellDict[depth][key].density

                    if adjacents[0] in otherCells:
                        if otherCells[adjacents[0]].density > threshold:
                            newCells[pstn] = remainHomoCells.pop(pstn)

                    elif adjacents[1] in otherCells:
                        if otherCells[adjacents[1]].density > threshold:
                            newCells[pstn] = remainHomoCells.pop(pstn)

                elif pstn in otherCells:
                    isFrozen = False
            
    # далее проверяются соседи Неймана, и подходящие добавляются в словарь
    # newCells; так же будет проверена текущая клетка, является ли она
    # замороженной (isFrozen) или остается активной 
            
            neumanNgbrs = self.activeCellDict[depth][key].getNeighbors()
            for pstn in neumanNgbrs:
                if pstn in remainHomoCells:
                    newCells[pstn] = remainHomoCells.pop(pstn)
                elif pstn in otherCells:
                    isFrozen = False
            
            # если подклеткам этой клетки далее некуда распространяться, то
            # эта клетка должна быть перемещена в словарь passiveCellDict:
            self.activeCellDict[depth].update(newCells)
            if isFrozen:
                self.passiveCellDict[depth][key] = self.activeCellDict[depth].pop(key)

def getAdjacentPair(firstPositin, secondPosition):
    """ функция получает в качестве аргумента позиции двух клеток, соседних 
    по гиперребру и возвращает позиции двух смежных с ними клеток """
    
    adjacent_1 = list(firstPositin)
    adjacent_2 = list(secondPosition)
    indices = list()
    for i in range(len(firstPositin)): 
        if firstPositin[i] != secondPosition[i]:
            indices.append(i)
    adjacent_1[indices[0]] = secondPosition[indices[0]]
    adjacent_2[indices[0]] = firstPositin[indices[1]]
    return tuple(adjacent_1), tuple(adjacent_2)

def expand(tempCellList, remainCells, cellList, otherCells, CONSTANTS):
    """ рекурсивная функция, расширяющая список клеток cellList кластера;
    в cellList содержатся клетки кластера, соседи которых уже проверены;
    tempCellList содержит клетки класстера, у которых могут быть соседи 
    в списке remainCells; в результате выполнения функции в списке cellList
    будут агрегированны все клетки кластера данной глубины, а из списка
    remainCells будут удалены все клетки, которые вошли в этот кластер """

    # условие выхода из рекурсии: исчерпываются клетки
    # у которых не проверены соседи
    if len(tempCellList) == 0: 
        return remainCells, cellList
    else:
        newTempCellList = []    # новый список который будет передан  
                                # рекурсивной функции вместо tempCellsList
        for cell in tempCellList:
            # перебираем клетки у которых ищем соседей Неймана
            for key in cell.getNeighbors(): # перебираем позиции соседей
                # если эта соседняя клетка находится в словаре remainCells:
                if key in remainCells.keys():
                    newTempCellList.append(remainCells.pop(key))
            # перебераем клетки, у которых ищем соседей по гиперребрам
            for key in cell.getExtraNeighbors():
                if key in remainCells.keys():
            # если такая клетка находится, то проверяем смежные с ними клетки
            # и объединяем в один кластер, если смежная клетка плотная
                    adjacents = getAdjacentPair(cell.getGlobalPosition(CONSTANTS), key)                    
                    threshold = CONSTANTS["ADJACENCY"] * cell.density

                    if adjacents[0] in otherCells:
                        if otherCells[adjacents[0]].density > threshold:
                            newTempCellList.append(remainCells.pop(key))

                    elif adjacents[1] in otherCells:
                        if otherCells[adjacents[1]].density > threshold:
                            newTempCellList.append(remainCells.pop(key))


            cellList.append(cell)   # список клеток данного кластера пополняется
        return expand(newTempCellList, remainCells, cellList, otherCells, CONSTANTS)

def showClusterCells(clusters, maxDepth, CONSTANTS, various=False):
    """Рисует рспределение класторов
    clusters: dict - клатсеры
    maxDepth: int - глубина самых мелких клеток (depth+2)
    """
    DENOMINATOR = CONSTANTS["DENOMINATOR"]
    def fillGrid(cell, grid, value):
        d = DENOMINATOR ** (maxDepth - cell.depth)
        gridPosition = cell.globalPosition[0] * d, cell.globalPosition[1] * d
        counter = 0
        for i in range(d):
            for j in range(d):
                counter += 1
                x = gridPosition[0]+i
                y = gridPosition[1]+j
                grid[x, y] = value
        return counter
    
    lenght = DENOMINATOR ** maxDepth
    grid = np.zeros((lenght, lenght))
    
    for depth in range(1, maxDepth):
        for clusterNumber in clusters:
            cluster = clusters[clusterNumber]
            value = clusterNumber if various else 1 
            if cluster.activeCellDict.get(depth, None):
                for pstn in cluster.activeCellDict.get(depth, None):
                    cell = cluster.activeCellDict[depth][pstn]
                    # print(f'active |{pstn}, depth: {cell.depth}')
                    fillGrid(cell, grid, value)

                for pstn in cluster.passiveCellDict.get(depth, None):
                    cell = cluster.passiveCellDict[depth][pstn]
                    # print(f'passive {pstn}, depth: {cell.depth}')
                    fillGrid(cell, grid, value)

    plt.xticks([])   # убираем
    plt.yticks([])   # деления
    fig = plt.gcf()
    fig.set_size_inches(8, 8)
    plt.imshow(grid.transpose(), origin='lower', cmap='gist_earth', alpha=0.8)


def clustering(X, CONSTANTS):
    N = len(X)
    # средняя плотность точек:
    meanDens =  N / CONSTANTS["SIZE"] ** (CONSTANTS["DIMENSION"] - 1)
    
    # создаем список с номерами всех экземпляров
    sampleList = [i for i in range(N)]
    # этот список будет также списком объектов клетки нулевой глубины

    progenitor = HeadCell(CONSTANTS)         # создается объект основной клетки
    progenitor.filling(sampleList, X, CONSTANTS)           # заполнение клетки
    # на этом этапе только создается и заполняется клетка нулевой глубины

    # Работа на первом уровне
    # ----------------------------------------------------------------------------

    # Получаем списки однородных и неоднородных непустых клеток первого уровня
    # в словарях `homoCells_1` и `otherCells_1`

    depth = 1                               # переходим на первый уровень
    # создаются и заполняются клетки, разбивающие основную клетку progenitor:
    # эти клетки [уровня 1] будут доступны через словарь otherCells_1
    otherCells_1 = progenitor.disperse(X, CONSTANTS) # (клетки первого уровня)

    # создаются и заполняются клетки, разбивающие клетки первого уровня:
    for key in otherCells_1:
        otherCells_1[key].disperse(X, CONSTANTS)
        # print(key, otherCells_1[key].density)

    # по статистическим данным данным, получаемым по клеткам второго уровня,
    # распределяются клетки первого уровня, словарь otherCells_1 обновляется:
    otherCells_1, homoCells_1 = progenitor.separation(depth, CONSTANTS)
    # otherCells_1 - прочие непустые клетки, homoCells_1 - однородные клетки

    # print(f'depth = 1\nЧисло однородных клеток:   {len(homoCells_1)}')
    # print(f'Число неоднородных клеток: {len(otherCells_1)}')


    # Работа на первом уровне. Образование кластеров из смежных клеток первого уровня.
    # Разделение клеток кластеров на активные (внешние) и пассивные (внутренние). 

    # класстеры будут содержаться в словаре clusters, и доступны по ключу номера кластера
    clusters = dict()
    clusterCounter = 0

    remainCells = copy.copy(homoCells_1)    # этот словарь будут исчерпываться в цикле

    # смежные однородные клетки первого уровня объединяются, образуя кластеры

    while len(remainCells) > 0:
        cellList = []       # здесь соберутся клетки одного кластера
        cell = remainCells.popitem()[1]        # извлекаем и удаляем клетку

        tempCellList = []   # список соберет соседние клетки клетки cell
        for key in cell.getNeighbors():     # перебираем список позиций соседей
            # если эта соседняя клетка находится в словаре remainCells, то
            # перемещаем ее из remainCells в список соседних клеток tempCellList
            if key in remainCells.keys():   
                tempCellList.append(remainCells.pop(key))
        cellList.append(cell) # cell становится первой клеткой новго кластера
        expand(tempCellList, remainCells, cellList, otherCells_1, CONSTANTS)
        # в результате в списке cellList будут собраны всек клетки, которые
        # находятся с клеткой cell в одном кластере
        clusterCounter += 1
        clusters[clusterCounter] = Cluster(cellList, depth, CONSTANTS)

    # print(f'Число кластеров: {len(clusters)}')

    # после того как кластеры образованы, нужно собрать их внешние (активные) клетки:
    for nmbr in clusters:
        clusters[nmbr].collectOuterSubCells(depth, X, CONSTANTS)
        # print(f'    active: {len(clusters[nmbr].activeCellDict[depth+1])}')


    # Работа на втором уровне
    # ----------------------------------------------------------------------------

    depth = 2       # переходим на второй уровень
    # собираем в словарь otherCells_2 клетки второго уровня из неоднородных
    # клеток первого уровня (по ключам, соответствующим глобальной позиции)
    otherCells_2 = dict()
    for key in otherCells_1:
        for k in otherCells_1[key].subCells:
            otherCells_2[otherCells_1[key].subCells[k].globalPosition] = otherCells_1[key].subCells[k]

    # создаются и заполняются клетки, разбивающие клетки второго уровня:
    for key in otherCells_2:
        otherCells_2[key].disperse(X, CONSTANTS)
    # по статистическим данным данным, получаемым по клеткам третьего уровня,
    # распределяются клетки второго уровня:
    for key in otherCells_1:
        otherCells_1[key].separation(depth, CONSTANTS)

    homoCells_2 = dict()    # будет содержать однородные клетки
    # словарь otherCells_2 обновляется, создается словарь с однородными клетками:
    for key in otherCells_1:
        for k in otherCells_1[key].subCells:
            otherCells_2[otherCells_1[key].subCells[k].globalPosition] = otherCells_1[key].subCells[k]
        for k in otherCells_1[key].homoCells:
            homoCells_2[otherCells_1[key].homoCells[k].globalPosition] = otherCells_1[key].homoCells[k]

    # print('otherCells_2', len(otherCells_2))
    # print('homoCells_2', len(homoCells_2))
    # print('clusters:', clusters)

    # showClusterCells(clusters, depth+2, CONSTANTS, various=False)


    # расширяем кластеры на однородные клетки второго уровня:
    while True:
        before = len(homoCells_2)
        for nmbr in clusters:
            clusters[nmbr].expand(homoCells_2, otherCells_2, depth, CONSTANTS)
        if before == len(homoCells_2):
            break

    # Работа на втором уровне. Из оставшихся клеток второго уровня выбирается 
    # произвольная клетка, которая образует новый кластер, после чего этот кластер 
    # тут же начинает расширяться на соседние однородные клетки второго уровня. 
    # Затем снова выбирается произвольная клетка из оставшихся и все повторяется 
    # до исчерпания всех оставшихся клеток

    remainCells = copy.copy(homoCells_2)

    # смежные однородные клетки второго уровня объединяются, образуя кластеры

    while len(remainCells) > 0:
        cellList = []       # здесь соберутся клетки одного кластера
        cell = remainCells.popitem()[1]        # извлекаем и удаляем клетку

        tempCellList = []   # список соберет соседние клетки клетки cell
        for key in cell.getNeighbors():     # перебираем список позиций соседей
            # если эта соседняя клетка находится в словаре remainCells, то
            # перемещаем ее из remainCells в список соседних клеток tempCellList
            if key in remainCells.keys():   
                tempCellList.append(remainCells.pop(key))
        cellList.append(cell) # cell становится первой клеткой новго кластера
        expand(tempCellList, remainCells, cellList, otherCells_2, CONSTANTS)
        # в результате в списке cellList будут собраны всек клетки, которые
        # находятся с клеткой cell в одном кластере
        clusterCounter += 1
        clusters[clusterCounter] = Cluster(cellList, depth, CONSTANTS)

    for cluster in clusters.values():
        cluster.collectOuterSubCells(depth, X, CONSTANTS)
        # print(f'    active: {len(cluster.activeCellDict[depth+1])}')

    # showClusterCells(clusters, depth+2, CONSTANTS, True)

    # Работа на втором уровне
    # ----------------------------------------------------------------------------

    depth = 3       # переходим на третий уровень
    # собираем в словарь otherCells_3 клетки третьего уровня из неоднородных
    # клеток второго уровня (по ключам, соответствующим глобальной позиции)
    otherCells_3 = dict()
    for key in otherCells_2:
        for k in otherCells_2[key].subCells:
            otherCells_3[otherCells_2[key].subCells[k].globalPosition] = otherCells_2[key].subCells[k]

    # создаются и заполняются клетки, разбивающие клетки третьего уровня:
    for key in otherCells_3:
        otherCells_3[key].disperse(X, CONSTANTS)
    # по статистическим данным данным, получаемым по клеткам третьего уровня,
    # распределяются клетки второго уровня:
    for cell in otherCells_2.values():
        cell.separation(depth, CONSTANTS)

    homoCells_3 = dict()
    # словарь otherCells_3 обновляется, создается словарь с однородными клетками:
    for key in otherCells_2:
        for k in otherCells_2[key].subCells:
            otherCells_3[otherCells_2[key].subCells[k].globalPosition] = otherCells_2[key].subCells[k]
        for k in otherCells_2[key].homoCells:
            homoCells_3[otherCells_2[key].homoCells[k].globalPosition] = otherCells_2[key].homoCells[k]

    # print('otherCells_3', len(otherCells_3))
    # print('homoCells_3', len(homoCells_3))
    # print('clusters:', clusters)

    # showClusterCells(clusters, depth+2, False, **CONSTANTS)

    # расширяем кластеры на однородные клетки третьего уровня:
    while True:
        before = len(homoCells_3)
        for nmbr in clusters:
            clusters[nmbr].expand(homoCells_3, otherCells_3, depth, CONSTANTS)
        if before == len(homoCells_3):
            break

    # print('homoCells_2', len(homoCells_3))
    # print('otherCells_2', len(otherCells_3))

    # showClusterCells(clusters, depth+2, True, **CONSTANTS)

    return clusters