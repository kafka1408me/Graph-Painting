import sys
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

from random import randint, shuffle
from PyQt5.QtCore import pyqtSignal, QObject, QThread
from PyQt5.QtGui import QIntValidator, QValidator
from PyQt5.QtWidgets import QApplication, QLabel, QMessageBox, QWidget, QMainWindow, QHBoxLayout, QVBoxLayout, \
    QLineEdit, QPushButton, QRadioButton, QTabWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

import graphs
import random
import graph_painting as MyGraph
import math
import enum

# --------------------------------------------------------------------

DEFAULT_NODE_COLOR = '#a7a7a7'

STOP_COND_MAX_GENERATION_STR = 'Выполнять до макс. поколения'
STOP_COND_PAINT_FOUND_STR    = 'Вернуть найденную раскраску сразу'

ALGRORITHM_GA_STR    = 'Генетический алгоритм'
ALGORITHM_IMMUNE_STR = 'Имунный алгоритм'

# problem constants:
HARD_CONSTRAINT_PENALTY = 10  # the penalty factor for a hard-constraint violation

class ProbablyValidator(QValidator):
    def validate(self, input_text, pos):
        try:
            val = float(input_text)
        except:
            return QValidator.Invalid, '',  pos

        if val > 0 and val < 1:
            return QValidator.Acceptable, input_text,  pos
        if val == 0:
            return QValidator.Intermediate, input_text,  pos
        return QValidator.Invalid, input_text,  pos

    pass

class MaxIntValidator(QValidator):

    def setMaxValue(self, max_value):
        self.max_value = max_value
        pass

    def validate(self, input_text, pos):
        if not input_text:
            return QValidator.Intermediate, input_text, pos
        try:
            val = int(input_text)
        except:
            return QValidator.Invalid, '', pos

        if val >= 0 and val < self.max_value:
            return QValidator.Acceptable, input_text,  pos
        return QValidator.Invalid, input_text, pos

    pass

def createPopulation(count_chromosomes, count_genes, max_colors):
    new_population = []
    for i in range(count_chromosomes):
        chromosome = []
        for j in range(count_genes):
            chromosome.append(random.randint(0, max_colors - 1))
            pass
        new_population.append(chromosome)
        pass
    return new_population

def mutation(chromosome, mut_percent, max_colors):
    if random.randint(0,100) < mut_percent:
        mut_indx = random.randint(0, len(chromosome) - 1)
        chromosome[mut_indx] = random.randint(0, max_colors - 1)
        pass
    pass

def createIndividual(count_genes, max_colors):
    max_color_indx = max_colors - 1
    return [random.randint(0, max_color_indx) for _ in range(count_genes)]

class StopCondition(enum.IntEnum):
    StopIfMaxGeneration = 0
    StopIfPaintFound    = 1

class FindPainting(QObject):

    foundPaintingResult = pyqtSignal(list)

    isStopped = False

    def __init__(self):
        super().__init__()

    def immuneAlgorithm(self, population_size, max_generations, max_colors, p_stay_percent, updated_percent, graph, stop_cond):
        gcp = graphs.GraphColoringProblem(graph=graph, hardConstraintPenalty=HARD_CONSTRAINT_PENALTY)
        count_genes = len(gcp)
        population = createPopulation(population_size, count_genes, max_colors)
        count_stay_p = int(p_stay_percent * population_size)
        count_updated = int(updated_percent * population_size)
        count_deleted_p = population_size - count_stay_p

        average_clones_for_one = population_size / count_stay_p
        p = 1.0 / count_genes
        population_indxs = [i for i in range(population_size)]
        list_indxs = [i for i in range(count_genes)]

        utilities = [(i, gcp.getCost(population[i])) for i in range(population_size)]
        utilities = sorted(utilities, key=lambda el: el[1])[:count_stay_p]

        for generation in range(max_generations):
            if self.isStopped:
                break
            print("Generation", generation)
            sum_utilities = 0
            for el in utilities:
                sum_utilities += el[1]

            new_population = []
            average_utility = sum_utilities / count_stay_p
            i = 0
            count_need_cloning = count_deleted_p
            utilities_2 = []
            while i < count_stay_p:
                count_clones = 0
                utility = utilities[i][1]
                if count_need_cloning > 0:
                    utility_relative = (average_utility / utility) * average_clones_for_one
                    count_clones = round(utility_relative)
                    count_need_cloning -= count_clones

                    if count_need_cloning < 0:
                        count_clones += count_need_cloning

                count_clones += 1
                utilities_2 += [(count_clones, utility)]

                new_population += [population[i] for _ in range(count_clones)]
                i += 1
                pass
           # print('new_population size =', len(new_population))

            # Мутация
            i = 0
            for count_antibodies, utility in utilities_2:
                mut = math.exp(p * utility)
                for k in range(count_antibodies):
                    count_changed_genes = math.ceil(mut)
                    if count_changed_genes > count_genes:
                        count_changed_genes = count_genes
                    random.shuffle(list_indxs)
                    for t in range(count_changed_genes):
                        population[i][list_indxs[t]] = random.randint(0, max_colors - 1)
                        pass
                    i += 1
                    pass
                pass
            population = new_population

            utilities = [(i, gcp.getCost(population[i])) for i in range(population_size)]
            utilities = sorted(utilities, key=lambda el: el[1])[:count_stay_p]

            if stop_cond == StopCondition.StopIfPaintFound:
                best = population[utilities[0][0]]
                if gcp.getViolationsCount(best) == 0:
                    return best

            if count_updated != 0:
                updated_indxs = random.sample(population_indxs, count_updated)
                for indx in updated_indxs:
                    population[indx] = createIndividual(count_genes, max_colors)

        best = population[utilities[0][0]]
        if gcp.getViolationsCount(best) == 0:
            return best
        return None



    def gaAlgorithm(self, population_size, max_generations, crossover_percent, mutation_percent, max_colors, graph, stop_cond):
        gcp = graphs.GraphColoringProblem(graph=graph, hardConstraintPenalty=HARD_CONSTRAINT_PENALTY)
        count_chromosome = len(gcp)
        population = createPopulation(population_size, count_chromosome, max_colors)

        def findBest(fitnesses=None):
            min_cost = math.inf
            best = None
            if fitnesses is None:
                fitnesses = [gcp.getCost(chromosome) for chromosome in population]

            for i in range(population_size):
                if fitnesses[i] < min_cost:
                    best = population[i]
                    pass
                pass

            if gcp.getViolationsCount(best) == 0:
                return best
            return None


        for generation in range(max_generations):
            if self.isStopped:
                break
            print('Generation', generation + 1)

            fitnesses = [gcp.getCost(chromosome) for chromosome in population]

            if stop_cond == StopCondition.StopIfPaintFound:
                best = findBest(fitnesses)
                if best is not None:
                    return best

            # турнирный отбор
            selected_indxs = []
            for j in range(population_size * 2):
                indx1 = random.randint(0, population_size - 1)
                indx2 = random.randint(0, population_size - 1)
                selected_indxs.append(
                    indx1 if fitnesses[indx1] < fitnesses[indx2] else indx2)
                pass
            # скрещивание
            new_population = []
            for j in range(0, population_size, 2):
                parent1_indx = selected_indxs[j]
                parent2_indx = selected_indxs[j + 1]
                child_1 = None
                child_2 = None
                if parent1_indx != parent2_indx and random.randint(0, 100) >= crossover_percent:
                    point_split = random.randint(1, count_chromosome - 1)
                    child_1 = population[parent1_indx][:point_split] + population[parent2_indx][point_split:]
                    child_2 = population[parent2_indx][:point_split] + population[parent1_indx][point_split:]
                    pass
                else:
                    child_1 = population[parent1_indx].copy()
                    child_2 = population[parent2_indx].copy()
                    pass
                new_population.append(child_1)
                new_population.append(child_2)
                pass
            population = new_population
            # mutation
            for chromosome in new_population:
                mutation(chromosome, mutation_percent, max_colors)
                pass
            pass

        min_cost = math.inf
        best = None
        for chromosome in population:
            cost = gcp.getCost(chromosome)
            if cost < min_cost:
                best = chromosome
                pass
            pass

        if gcp.getViolationsCount(best) == 0:
            return best
        return None

    def interrupt(self):
        self.isStopped = True
        pass

    def findGraphPaintingGA(self, population_size, max_generations, crossover_percent, mutation_percent, max_colors, graph, stop_cond):
        self.isStopped = False
        best = self.gaAlgorithm(population_size=population_size, max_generations=max_generations,
                                   crossover_percent=crossover_percent,
                                   mutation_percent=mutation_percent,
                                   max_colors=max_colors, graph=graph, stop_cond=stop_cond)
        if best is None:
            self.foundPaintingResult.emit([])
        else:
            self.foundPaintingResult.emit(best)
        pass

    def findPaintingImmune(self, population_size, max_generations, max_colors, p_stay_percent, updated_percent, graph, stop_cond):
        self.isStopped = False
        best = self.immuneAlgorithm(population_size=population_size, max_generations=max_generations,
                                    max_colors=max_colors, p_stay_percent=p_stay_percent,
                                    updated_percent=updated_percent, graph=graph, stop_cond=stop_cond)

        if best is None:
            self.foundPaintingResult.emit([])
        else:
            self.foundPaintingResult.emit(best)
        pass
    pass

class WidgetGA(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        self.mutation_probably_input = Input()
        self.mutation_probably_input.setTitle("Вероятность мутации")
        self.mutation_probably_input.setText("0.1")
        self.mutation_probably_input.setValidator(ProbablyValidator())

        self.crossover_probably_input = Input()
        self.crossover_probably_input.setTitle("Вероятность кроссовера")
        self.crossover_probably_input.setText("0.9")
        self.crossover_probably_input.setValidator(ProbablyValidator())

        layout.addWidget(self.mutation_probably_input)
        layout.addWidget(self.crossover_probably_input)
        self.setLayout(layout)

    def getMutation(self):
        return float(self.mutation_probably_input.text())

    def getCrossover(self):
        return float(self.crossover_probably_input.text())

class WidgetImmune(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()
        self.best_part_antibodies_input  = Input()
        self.best_part_antibodies_input.setTitle("Выбираемые лучшие антитела")
        self.best_part_antibodies_input.setText("0.3")
        self.best_part_antibodies_input.setValidator(ProbablyValidator())

        self.updated_antibodies_input = Input()
        self.updated_antibodies_input .setTitle("Обновляемые антитела")
        self.updated_antibodies_input .setText("0.1")
        self.updated_antibodies_input .setValidator(ProbablyValidator())

        layout.addWidget(self.best_part_antibodies_input)
        layout.addWidget(self.updated_antibodies_input)
        self.setLayout(layout)

    def getBestAntibodiesPart(self):
        return float(self.best_part_antibodies_input.text())

    def getUpdatedPart(self):
        return float(self.updated_antibodies_input.text())


class MainWidget(QWidget):

    find_painting = FindPainting()

    sendParamsForGA = pyqtSignal(int, int, float, float, int, nx.Graph, StopCondition)

    sendParamsForImmune = pyqtSignal(int, int, int, float, float, nx.Graph, StopCondition)

    def __init__(self):
        super().__init__()
        self.initUI()

        self.thread_finder_painting = QThread()
        self.find_painting.moveToThread(self.thread_finder_painting)
        self.thread_finder_painting.start()

        self.find_painting.foundPaintingResult.connect(self.paintingResult)
        self.sendParamsForGA.connect(self.find_painting.findGraphPaintingGA)
        self.sendParamsForImmune.connect(self.find_painting.findPaintingImmune)
        pass

    def initUI(self):
        main_layout = QHBoxLayout()
        sidebar_layout = QVBoxLayout()
        sidebar_layout.setSpacing(5)

        self.canvas = MplCanvas(self, 5, 4)

        self.graph_size_input = Input()
        self.graph_size_input.setTitle("Количество вершин")
        self.graph_size_input.setText("12")
        self.graph_size_input.setValidator(QIntValidator())

        self.graph_size_input.textChanged.connect(self.setMaxEdgeCount)

        self.edge_validator = MaxIntValidator()
        self.setMaxEdgeCount()

        self.number_of_edges_input = Input()
        self.number_of_edges_input.setTitle("Количество рёбер")
        self.number_of_edges_input.setText("17")
        self.number_of_edges_input.setValidator(self.edge_validator)

        self.population_size = Input()
        self.population_size.setTitle("Размер популяции")
        self.population_size.setText("100")
        self.population_size.setValidator(QIntValidator())

        self.max_generations = Input()
        self.max_generations.setTitle("Количество поколений")
        self.max_generations.setText("100")
        self.max_generations.setValidator(QIntValidator())

        self.max_count_colors_input = Input()
        self.max_count_colors_input.setTitle("Максимальное количество цветов")
        self.max_count_colors_input.setText("5")
        self.max_count_colors_input.setValidator(QIntValidator())

        self.generate_button = QPushButton('Сгенерировать', self)
        self.generate_button.clicked.connect(self.generate)

        self.run_button = QPushButton('Найти раскраску', self)
        self.run_button.clicked.connect(self.run)

        self.result_label = QLabel()
        self.result_label.setWordWrap(True)

        self.interrupt_btn = QPushButton('Прервать вычисление', self)
        self.interrupt_btn.clicked.connect(self.interruptFindingPainting)
        self.interrupt_btn.setVisible(False)


        main_layout.addWidget(self.canvas)
        main_layout.setStretch(0, 10)
        main_layout.addItem(sidebar_layout)

        self.choice_algorithm_stop = ChoiceBetween([STOP_COND_MAX_GENERATION_STR, STOP_COND_PAINT_FOUND_STR])

        self.algorithms_tab = QTabWidget()
        self.algorithms_tab.setMinimumWidth(240)
        self.ga_widget = WidgetGA()
        self.immune_widget = WidgetImmune()
        self.algorithms_tab.addTab(self.ga_widget, ALGRORITHM_GA_STR)
        self.algorithms_tab.addTab(self.immune_widget, ALGORITHM_IMMUNE_STR)

        self.algorithms_tab.setStyleSheet("""QTabWidget::pane { /* The tab widget frame */
    border-top: 2px solid #C2C7CB;
    position: absolute;
    top: -0.5em;
    border-bottom: 2px solid #C2C7CB;
    background-color: #dfdfdf;
}

QTabWidget::tab-bar {
    alignment: left;
}


/* Style the tab using the tab sub-control. Note that
    it reads QTabBar _not_ QTabWidget */
QTabBar::tab {
    background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                                stop: 0 #E1E1E1, stop: 0.4 #DDDDDD,
                                stop: 0.5 #D8D8D8, stop: 1.0 #D3D3D3);
    border: 2px solid #C4C4C3;
    border-bottom-color: #C2C7CB; /* same as the pane color */
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    min-width: 8ex;
    padding: 2px;
}

QTabBar::tab:selected, QTabBar::tab:hover {
    background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                                stop: 0 #fafafa, stop: 0.4 #f4f4f4,
                                stop: 0.5 #e7e7e7, stop: 1.0 #fafafa);
}

QTabBar::tab:selected {
    border-color: #9B9B9B;
    border-bottom-color: #C2C7CB; /* same as pane color */
}""")

        sidebar_layout.addWidget(self.algorithms_tab)
        sidebar_layout.addWidget(self.graph_size_input)
        sidebar_layout.addWidget(self.number_of_edges_input)
        sidebar_layout.addWidget(self.population_size)
        sidebar_layout.addWidget(self.max_generations)
        sidebar_layout.addWidget(self.max_count_colors_input)
        sidebar_layout.addWidget(self.choice_algorithm_stop)
        sidebar_layout.addWidget(self.generate_button)
        sidebar_layout.addWidget(self.run_button)
        sidebar_layout.addWidget(self.result_label)
        sidebar_layout.addWidget(self.interrupt_btn)
        sidebar_layout.addStretch(10)

        self.setLayout(main_layout)
        pass

    def paintingResult(self, res):
        if len(res) != 0:
            self.draw(res)
        else:
            self.result_label.setText("<font color='red'>Не удалось найти раскраску</font>")
            self.draw()
        self.setCalculateState(False)
        pass

    def setCalculateState(self, isCalculating: bool) -> None:
        isSetGraphParamsEnabled = not isCalculating
        self.graph_size_input.setEnabled(isSetGraphParamsEnabled)
        self.number_of_edges_input.setEnabled(isSetGraphParamsEnabled)
        self.algorithms_tab.setEnabled(isSetGraphParamsEnabled)
        self.max_count_colors_input.setEnabled(isSetGraphParamsEnabled)
        self.population_size.setEnabled(isSetGraphParamsEnabled)
        self.max_generations.setEnabled(isSetGraphParamsEnabled)
        self.generate_button.setEnabled(isSetGraphParamsEnabled)
        self.run_button.setEnabled(isSetGraphParamsEnabled)
        self.choice_algorithm_stop.setEnabled(isSetGraphParamsEnabled)
        self.interrupt_btn.setVisible(isCalculating)
        pass

    def interruptFindingPainting(self):
        self.find_painting.interrupt()

    def run(self):

        try:
            max_colors = int(self.max_count_colors_input.text())
            if not max_colors:
                raise Exception()
        except:
            self.showError('Введите максимальное количество цветов для раскраски')
            return

        try:
            population_size = int(self.population_size.text())
            if not population_size:
                raise Exception()
        except:
            self.showError('Введите размер популяции')
            return

        try:
            max_generations = int(self.max_generations.text())
            if not max_generations:
                raise Exception()
        except:
            self.showError('Введите количество поколений')
            return


        algorithm_stop_cond = self.choice_algorithm_stop.getCheckedStr()
        if algorithm_stop_cond == STOP_COND_MAX_GENERATION_STR:
            algorithm_stop_cond = StopCondition.StopIfMaxGeneration
        else:
            algorithm_stop_cond = StopCondition.StopIfPaintFound


        if self.algorithms_tab.currentIndex() == 0: # GA
            try:
                mutation_probability = self.ga_widget.getMutation()
            except:
                self.showError('Введите вероятность мутации')
                return

            try:
                crossover_probability = self.ga_widget.getCrossover()
            except:
                self.showError('Введите вероятность кроссовера')
                return
            self.sendParamsForGA.emit(population_size, max_generations,
                                      crossover_probability * 100,
                                      mutation_probability * 100,
                                      max_colors, self.graph.copy(), algorithm_stop_cond)
            algorithm_name = 'генетическим алгоритмом'
        else: # Immune
            try:
                best_part_antibodies = self.immune_widget.getBestAntibodiesPart()
            except:
                self.showError('Введите вероятность мутации')
                return
            try:
                updated_part = self.immune_widget.getUpdatedPart()
            except:
                self.showError('Введите вероятность кроссовера')
                return
            self.sendParamsForImmune.emit(population_size, max_generations,
                                     max_colors, best_part_antibodies, updated_part, self.graph, algorithm_stop_cond)
            algorithm_name = 'имунным алгоритмом'

        self.result_label.setText("<font color='#0c7af7'>Идет поиск раскраски {}..</font>".format(algorithm_name))
        self.setCalculateState(True)
        pass

    def setMaxEdgeCount(self):
        try:
            count_nodes = int(self.graph_size_input.text())
        except:
            return

        max_edges = count_nodes*(count_nodes-1) // 2
        self.edge_validator.setMaxValue(max_edges)

        try:
            current_count_edge = int(self.number_of_edges_input.text())
        except:
            return

        if current_count_edge > max_edges:
            self.number_of_edges_input.setText('0')

        pass


    def generate(self):
        self.result_label.setText("")

        try:
            count_nodes = int(self.graph_size_input.text())
        except ValueError:
            self.showError("Введите размер графа")
            return

        try:
            count_edges = int(self.number_of_edges_input.text())
        except ValueError:
            self.showError("Введите количество рёбер")
            return

        k = count_nodes - 1
        num_node = 0
        list_adjacency = []
        while k != 0:
            list_adjacency += [x for x in range(num_node * count_nodes + num_node + 1, count_nodes * (num_node + 1))]
            #    print(num_node*(count_nodes + 1) + num_node + 1, count_nodes*(num_node + 1))
            num_node += 1
            k -= 1
            pass

        random.shuffle(list_adjacency)

        list_edges = []
        for x in list_adjacency[:count_edges]:
            list_edges += [(x // count_nodes, x % count_nodes)]
            pass

        list_other_nodes = [x // count_nodes for x in list_adjacency[count_edges:]]

        self.graph = nx.from_edgelist(list_edges)
        self.graph.add_nodes_from(list_other_nodes)

        self.draw()
        pass

    def draw(self, colorArrangement=None):
        self.canvas.clear()

        if colorArrangement is not None:

            # create a list of the unique colors in the arrangement:
            colorList = list(set(colorArrangement))

            self.result_label.setText("<font color='green'>Найдена раскраска с {} цветами</font>".format(len(colorList)))

            # create the actual colors for the integers in the color list:
            colors = plt.cm.rainbow(np.linspace(0, 1, len(colorList)))

            # iterate over the nodes, and give each one of them its corresponding color:
            colorMap = []
            for i in range(nx.number_of_nodes(self.graph)):
                colorMap.append(colors[colorList.index(colorArrangement[i])])
                pass
            pass
        else:
            colorMap = [DEFAULT_NODE_COLOR for _ in range(self.graph.number_of_nodes())]
            pass

        nx.draw_circular(self.graph, ax=self.canvas.axes, node_color=colorMap, with_labels=True)

        self.canvas.draw()
        pass

    def showError(self, text: str):
        error_dialog = QMessageBox()
        error_dialog.setText(text)
        error_dialog.adjustSize()
        error_dialog.exec()
        pass

    pass


class Input(QWidget):
    textChanged = pyqtSignal(str)

    def __init__(self):
        QWidget.__init__(self)
        self.initUI()

    def slot_textChanged(self, input_text):
        self.textChanged.emit()
        pass

    def initUI(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        self.title = QLabel()
        self.le = QLineEdit()

        self.le.textChanged.connect(self.textChanged)

        layout.addWidget(self.title)
        layout.addWidget(self.le)

    def text(self) -> str:
        return self.le.text()

    def setText(self, text: str):
        self.le.setText(text)

    def setTitle(self, text: str):
        self.title.setText(text)

    def setValidator(self, validator: QValidator):
        self.le.setValidator(validator)

    def setEnabled(self, enabled: bool) -> None:
        self.le.setEnabled(enabled)

    pass

class ChoiceBetween(QWidget):

    def __init__(self, variants):
        QWidget.__init__(self)
        self.initUI(variants)

    def initUI(self, variants):
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        btn = QRadioButton(variants[0])
        btn.setChecked(True)
        self.layout.addWidget(btn)
        for variant in variants[1:]:
            btn = QRadioButton(variant)
            self.layout.addWidget(btn)

    def getCheckedStr(self) -> str:
        for i in range(self.layout.count()):
            item = self.layout.itemAt(i).widget()
            if item.isChecked():
                return item.text()
        return ''

    def setEnabled(self, enabled: bool) -> None:
        for i in range(self.layout.count()):
            item = self.layout.itemAt(i).widget()
            item.setEnabled(enabled)



class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        super(MplCanvas, self).__init__(self.fig)
        pass

    def clear(self):
        self.axes.cla()
        pass

    pass


class Window(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setGeometry(500, 300, 1000, 700)
        self.setWindowTitle('Нахождение раскраски графа')

        mainWidget = MainWidget()
        mainWidget.generate()

        self.setCentralWidget(mainWidget)
        self.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Window()
    sys.exit(app.exec())
