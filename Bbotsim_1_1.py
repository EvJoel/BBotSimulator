import random as rd
import re
from math import *
import time
import math 
import matplotlib.pyplot as plt
from Tkinter import *
import Tkinter
import ttk 
import tkMessageBox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
#Global variables

test = [[[0, 0], [0.08533058638602291, -0.08434345170904825]], [[0.08533058638602291, -0.08434345170904825], [0.1248298766273177, 0.02571312418400408]], [[0.1248298766273177, 0.02571312418400408], [0.12482987662731772, 0.025713124184004107]], [[0.12482987662731772, 0.025713124184004107], [0.18872532857212487, 0.018942627380856808]], [[0.18872532857212487, 0.018942627380856808], [0.3932960503175464, -0.17052878609697728]]]


NUM_SEQ = 5
NUM_CHROMES = 5
NUM_GEN = 10
CURR_GEN = 10

CHROM_MUTATE_RATE = .2
GENE_MUTATE_RATE = .33

host = "10.201.109.60"
port = 1337

#controls the weights for the fitness function
k = [0.75, 0.15, 0.5, 0.5]


START = [0,0] #Markers for the starting area
GOAL = [0,1] #Markers for the finish area

###########################################################################################################################################################


def createGeneration():
    generation = [[[rd.randint(0,180), rd.randint(0,180), rd.randint(1,1)] for i in range(NUM_CHROMES)] for i in range(NUM_SEQ)] #Creates a generation with chromosomes ranging between (0,180) and (1,4)
    print "Created Generation with " + str(NUM_SEQ) + " sequences consisting of " + str(NUM_CHROMES) + " Chromosomes" 
    return generation
    
def printSequence(seq):
    for i in range(len(seq)): # prints a sequence in a readable way
        print seq[i]
        
def getDistance(x1,y1,x2,y2):
    distance = sqrt(((x2-x1)**2) + ((y2-y1)**2))
    return distance 

def mutateSequence(seq):
    global MUTATE_NUM
    new_seq = []
    for i in range(len(seq)):
        new_seq.append(list(seq[i]))

    for i in range(len(new_seq)):
        mutate_chance_chromo = 100 * CHROM_MUTATE_RATE #Numbers used to decide whether a chromosome in the sequence is mutated
        mutate_number_chromo = rd.randint(1,100)

        if(mutate_number_chromo <= mutate_chance_chromo):
            #if one is greater, it has the possibility to mutate
            for j in range(3):
                #loops through each gene
                pos_neg = rd.randint(0,100) #add or subtract from gene
                mutate_number_gene = rd.randint(0,100) #Numbers used to decide whether a gene in the sequence is mutated
                mutate_chance_gene = 100 * GENE_MUTATE_RATE
                MUTATE_NUM = rd.randint(1,100)  #amount gene is changed by
                if(mutate_number_gene <= mutate_chance_gene):
                    change = 0
                    if(j < 2):
                        if(mutate_chance_gene <= 50):
                            change =  MUTATE_NUM*.1
                        elif(mutate_chance_gene <= 70 and mutate_rate > 50):
                            change = MUTATE_NUM *.15
                        elif(mutate_chance_gene <= 80 and mutate_rate > 70):
                            change =  MUTATE_NUM*.20
                        elif(mutate_chance_gene <= 85 and mutate_rate > 80):
                            change =  MUTATE_NUM*.25
                        elif(mutate_chance_gene <= 90 and mutate_rate > 85):
                            change =  MUTATE_NUM*.30
                        elif(mutate_chance_gene <= 94 and mutate_rate > 90):
                            change =  MUTATE_NUM*.4   
                        elif(mutate_chance_gene <= 97 and mutate_rate > 94):
                            change =  MUTATE_NUM*.5
                        elif(mutate_chance_gene <= 99 and mutate_rate > 97):
                            change =  MUTATE_NUM*.6
                        else:
                            change =  MUTATE_NUM*.7

                        if(pos_neg >= 50):
                            #print new_seq[i][j]
                            new_seq[i][j]  = new_seq[i][j] + int(change)
                            #print "pos " +str(change)
                        else:
                            new_seq[i][j] = new_seq[i][j] - int(change)
                            #print "neg " +  str(change)

                        if(new_seq[i][j] > 180):
                            new_seq[i][j] = 180
                        if(new_seq[i][j] < 0):
                            new_seq[i][j] = 0          

                    else:
                        #used for mutating time genes
                        if(new_seq[i][j] ==1):
                            if(mutate_chance_gene <= 60):
                                new_seq[i][j] = 2
                            elif(mutate_chance_gene > 60 and mutate_rate <= 90):
                                new_seq[i][j] =  3
                            else:
                                new_seq[i][j] = 4

                        elif(new_seq[i][j] ==2):
                            if(mutate_chance_gene <= 40):
                                new_seq[i][j] = 1
                            elif(mutate_chance_gene > 40 and mutate_rate <= 80):
                                new_seq[i][j] =  3
                            else:
                                new_seq[i][j] = 4

                        elif(new_seq[i][j] ==3):

                            if(mutate_chance_gene <= 60):
                                new_seq[i][j] = 3
                            elif(mutate_chance_gene > 60 and mutate_rate <= 90):
                                new_seq[i][j] =  4
                            else:
                                new_seq[i][j] = 1  

                        else:
                            if(mutate_chance_gene <= 60):
                                new_seq[i][j] = 3
                            elif(mutate_chance_gene > 60 and mutate_rate <= 90):
                                new_seq[i][j] =  2
                            else:
                                new_seq[i][j] = 1                           


    return new_seq

def getFitness(dist_from_goal, dist_traveled, elap_time, dist_goal_final):
    fitness = (float(dist_from_goal) * k[0]) + (float(dist_traveled) * k[1]) + ((float(elap_time) * k[2])/10) + (float(dist_goal_final) * k[3])
    return fitness


def readFitness(gen_num):
    generation_new = [[0,0,0] for i in range(NUM_CHROMES)]
    current_gene = 0
    most_fit = 100000000000
    most_fit_num = 0
    
    import math
    import matplotlib.pyplot as plt
    
    
def calculuteBobotMotion(x0,y0,theta0,m1,m2,T,N): 
    # INPUT
    #     x0,y0  : (float) Coordinates of the center of the bobot 
    #     theta0 : (float) Direction Bobot is facing
    #     m1,m2  : (int) Servo values from 0 to 180 for left and right servo
    #     T      : (float) Number of seconds from 1 to 4 to run servos
    #     N      : (int) The number of coordinates you want to break the trajectory into
    #
    # OUTPUT
    #     POS      : ((float,float)) List of (x,y) coordinates from the start position to end position
    #     theta    : (float) Updated bobot direction
    #     distance : (float) distance in meters bobot traveled
    # 
    
    # CONSTANTS 
    # Units are in meters and seconds
    chassis_width = .12     
    wheel_diameter = .067
    
    # Calculuate distances the wheels turn from m1, m2
    # First determines the rate of rotation (rot per second), then distance = rotation*(2*pi*r)*T
    d1 = ((abs(m1)-90.0)/90.0)*(5.0/6.0) * math.pi*wheel_diameter * T
    d2 = ((abs(m2)-90.0)/90.0)*(5.0/6.0) * math.pi*wheel_diameter * T
        
    # Calculate the updated position
    POS = [[x0,y0]]
    if (d1 == d2):
        # Forward or backward case   
        delta_d = d1/N 
        for step in range(N):
            d = (step+1)*delta_d
            x = x0+d*math.cos(theta0)
            y = y0+d*math.sin(theta0)
            POS.append([x,y])
        theta = theta0
        distance = d1
        #print distance
    elif (d1>d2):
        # Right turn (forward and backward)
        theta = (d1-d2)/chassis_width
        delta_theta = theta / N
        radius = d1 + chassis_width/2.0       
        distance = theta*radius
        for step in range(1,N+1):
            d = (step)*delta_theta 
            x = x0 + math.cos(theta0-math.pi/2.0)*radius + math.cos(theta0 + math.pi/2.0 - step*delta_theta)*radius 
            y = y0 + math.sin(theta0-math.pi/2.0)*radius + math.sin(theta0 + math.pi/2.0 - step*delta_theta)*radius
            POS.append([x,y])
        distance = radius * theta 
            
    else: # d1<d2
        # Left turn (forward and backward)
        theta = (d2-d1)/chassis_width
        delta_theta = theta / N
        radius = d2 + chassis_width/2.0
        distance = theta * radius
        for step in range(N):
            d = step*delta_theta 
            x = x0 + math.cos(theta0+math.pi/2.0)*radius + math.cos(theta0 - math.pi/2.0 + step*delta_theta)*radius 
            y = y0 + math.sin(theta0+math.pi/2.0)*radius + math.sin(theta0 - math.pi/2.0 + step*delta_theta)*radius
            POS.append([x,y])      
        distance = radius * theta 
        
    return POS,theta,distance   

def checkDistFromGoal(position,dist_from_goal):
    for i in range(len(position)):
        if getDistance(position[i][0],position[i][1],GOAL[0],GOAL[1]) >= dist_from_goal:
            dist_from_goal = getDistance(position[i][0],position[i][1],GOAL[0],GOAL[1])
    return dist_from_goal
            
    
def createGraph(pos):
    plt.clf()
    #
    #plt.set_animated(True)
    print len(pos)
    x = []
    y = []
    for i in range(len(pos)):
        for j in range(len(pos[i])):
            print pos[i][j]
            x.append(pos[i][j][0])
            y.append(pos[i][j][1])
    
    plt.plot(x,y,'ro')
    plt.ylabel('Y')
    plt.xlabel('X')
    plt.axis([-1,1,-1,1])
    plt.axis('equal')
    plt.show()
    

    
    
class Simulator:
    def __init__(self,master,enableGraph):
        self.enableGraph = enableGraph
        self.fig = Figure(figsize=(5,4),dpi=100)
        self.graph = self.fig.add_subplot(111)        
        
        self.master = master
        self.frame = Frame(self.master)
        self.frame.pack()
    
    
        self.leftFrame = Frame(root)
        self.leftFrame.pack(side=LEFT,expand=1)
        
        self.rightFrame = Frame(root)
        self.rightFrame.pack(side = RIGHT)
        
        self.bottomFrame = Frame(self.leftFrame)
        self.bottomFrame.pack( side = BOTTOM )  
        
        self.bbFrame = Frame(self.bottomFrame)
        self.bbFrame.pack(side=BOTTOM)
        
        self.numGenLabel = Label(self.leftFrame, text="Number of Generations",)
        self.numGenLabel.pack(side=LEFT,fill=X)
        
        self.numGen = Entry(self.leftFrame, bd = 5)
        self.numGen.pack(side=RIGHT,fill=X)
        
        self.startButton = Button(self.bottomFrame, text = "Run Simulator", command = self.runMain)
        self.startButton.pack(side = TOP,fill=X)
        
        if enableGraph:
            self.graphFrame = FigureCanvasTkAgg(self.fig,self.rightFrame)
            self.graphFrame.show()
            self.graphFrame.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1) 
        else:
            self.graph = Canvas(self.rightFrame,width=400,height=400)
            self.graph.config(bg="BLUE")
            self.graph.create_text(200,200,text="This is a placeholder",fill="WHITE")
            self.graph.pack(side=RIGHT)            
        
        self.currentGen = Label(self.bottomFrame,text="Current Generation: ")
        self.currentGen.pack(side=LEFT,)
        
        self.genNum = Label(self.bottomFrame,text="0")
        self.genNum.pack(side=RIGHT,)
        
        
        self.progressBar = ttk.Progressbar(self.bbFrame,orient=HORIZONTAL, length=200, mode='determinate')
        self.progressBar.pack()
        
        #self.updateGUI()
    def runMain(self):
        num = self.numGen.get()
        if num:
            #self.updateGraph()
            self.runMainProgram(int(num))
        else:
            tkMessageBox.showinfo("Error", "Please enter the number of generations you wish to run") 
        
        
    def runMainProgram(self,num_gen):
        generation = createGeneration()
        dist_from_goal = 0
        current_seq = 0
        current_gen = 0
        new_gen = 0
        theta = 0
        gen_fitness = []
        distance = 0
        dist_total = 0
        POS = []
        x0 = 0
        y0 = 0
        total_time = 0
        original = generation
        first_gen_fitness = []
        fitness_scores = []
        position_total = []
        self.progressBar.config(maximum=num_gen*5)
        #time.sleep(5)
        compare_gen = num_gen
        start = time.time()
        
        while num_gen > new_gen:
            
            for i in range(NUM_SEQ):
                for j in range(NUM_CHROMES):
                    #print generation[i][j]
                    POS , theta, distance = calculuteBobotMotion(x0, y0, theta, generation[i][j][0], generation[i][j][1], generation[i][j][2], 100)
                    x0 = POS[len(POS)-1][0]
                    y0 = POS[len(POS)-1][1]
                    #print theta
                    #print x0
                    #print y0
                    position_total.append(POS)
                    #print position_total
                    dist_total += distance
                    distance = 0
                    dist_from_goal = checkDistFromGoal(POS, dist_from_goal)
                    total_time += generation[i][j][2]
                    #print total_time
                goal_dist = getDistance(x0, y0, GOAL[0], GOAL[1]) 
                gen_fitness.append([current_seq, getFitness(dist_from_goal, dist_total, total_time, goal_dist)])
                self.progressBar.step(1)
                self.progressBar.update()
                #print  "Sequence fitness: " + str(getFitness(dist_from_goal, dist_total, total_time, goal_dist))  
                #createGraph(position_total)
                if self.enableGraph:
                    self.getNewGraph(position_total)
                position_total = []
                x0 = 0
                y0 = 0
                theta = 0
                dist_total = 0
                total_time = 0 
                goal_dist = 0
                POS = []
                dist_from_goal = 0
                current_seq += 1
                #print self.graphFig.axes
                #time.sleep(2)
    
            most_fit = [0,100000000000]

            #print most_fit
            for i in range(len(gen_fitness)):
                if(gen_fitness[i][1] < most_fit[1]):
                    most_fit = gen_fitness[i]
                    
            mutated_gen = []
            #print most_fit
            #time.sleep(2)
            fitness_scores.append([(compare_gen-num_gen), most_fit[1]])
            
            #print generation
            #print gen_fitness
            
            #print "Most fit sequence: " + str(most_fit)
           # print mutated_gen
            #print most_fit
            if num_gen == compare_gen:
                first_gen_fitness = most_fit
            for i in range(NUM_SEQ):
                #print mutated_gen.append(mutateSequence(generation[most_fit[0]]))
                mutated_gen.append(mutateSequence(generation[most_fit[0]]))
            #print generation 
            generation = mutated_gen   
            #print generation
            num_gen -= 1
            current_seq = 0
            gen_fitness = []
            self.genNum.config(text=str(compare_gen-num_gen))
            self.genNum.update()
        print fitness_scores
        self.graphFitness(fitness_scores)
        
            
       
    def getNewGraph(self,pos):
        #plt.ion() 
        x = []
        y = []
        for i in range(len(pos)):
            for j in range(len(pos[i])):
                #print pos[i][j]
                x.append(pos[i][j][0])
                y.append(pos[i][j][1])
                
        self.graph.plot(x,y,'ro',)
        #graph.ylabel('Y')
        #graph.xlabel('X')
        self.graph.axis([-1,1,-1,1])
        self.graph.axis('equal')
        self.graphFrame.draw()
        self.graph.clear()
        
    def graphFitness(self,fitness_scores):
        x = []
        y = []
        for i in range(len(fitness_scores)):
            print fitness_scores[i]
            x.append(fitness_scores[i][0])
            y.append(fitness_scores[i][1])
                
        plt.plot(x,y)
        plt.ylabel('Fitness Score')
        plt.xlabel('Generation')
        plt.axis([0,len(fitness_scores),0,2],)
        #plt.axis('equal')
        plt.show()        
          
 
            
    
    
        
        
            

    
#####################################################################################################################################################################################
    
#Main execution loop 
root = Tk()

app = Simulator(root,False)


root.mainloop()


        
                
            
            
    

    