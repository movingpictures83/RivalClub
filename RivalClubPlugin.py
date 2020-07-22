from AffinityPropagation.AffinityPropagationPlugin import *
import math

class RivalClubPlugin(AffinityPropagationPlugin):
   def input(self, inputfile):
      filestuff = open(inputfile, 'r')
      # Populate two ADJ matrices
      # One with the original network
      # The other the same network, but all negatives become positives
      firstline = filestuff.readline()
      self.bacteriaORIG = firstline.split(',')
      if (self.bacteriaORIG.count('\"\"') != 0):
         self.bacteriaORIG.remove('\"\"')
      self.n = len(self.bacteriaORIG)
      self.ADJNONEG = []
      self.ADJORIG = []
      i = 0
      for line in filestuff:
         contents = line.split(',')
         self.ADJNONEG.append([])
         self.ADJORIG.append([])
         for j in range(self.n):
            value = float(contents[j+1])
            if (i != j and value != 0):
               self.ADJNONEG[i].append(abs(value))
               self.ADJORIG[i].append(value)
            else:
               self.ADJNONEG[i].append(0)
               self.ADJORIG[i].append(0)
         i += 1

   def doCluster(self, ADJMAT, BAC):
      self.ADJ = ADJMAT.copy()
      self.bacteria = BAC.copy()
      eps = 1e-8
      ap = AffinityPropagation(preference=0,damping=0.75,affinity='precomputed',convergence_iter=200, verbose=False)
      self.removeSingletonsAndPureVillains()
      #print("AFTER REMOVING SINGLETONS: "+str(self.ADJ))
      #print(self.bacteria)
      af = ap.fit(self.ADJ)
      self.cluster_centers_indices = af.cluster_centers_indices_
      self.labels = af.labels_
      self.n_clusters_ = len(self.cluster_centers_indices)
      #print(self.cluster_centers_indices)
      #print(self.labels)
      #y = input()
      clusters = []
      #pr`int(self.labels)
      #print(self.labels[i])
      # Mark singletons, will remove when printing
      cluster = 1
      for i in range(0, len(self.cluster_centers_indices)):
         cluster_size = 0
         thecluster = []
         for j in range(0, len(self.labels)):
            label = self.labels[j]
            if (label == cluster):
               cluster_size += 1
               thecluster.append(self.bacteria[j])
         if (cluster_size <= 1):
            self.cluster_centers_indices[i] = -1 # Mark
         else:
            clusters.append(thecluster)
         cluster += 1
      #print("CLUSTERS:"+str(clusters))
      return clusters
            

   def run(self):
      self.rivalclubs = [] 
      # CLUSTER NONEG NETWORK FIRST
      superclusters = self.doCluster(self.ADJNONEG, self.bacteriaORIG)
      print("TOTAL SUPERCLUSTERS: "+str(len(superclusters))) 
      #self.ADJNEW = self.ADJ.copy()
      #print(self.ADJNEW)
      #print("***")
      ## NOW CLUSTER THE SUPERCLUSTERS
      # USE ORIGINAL MATRIX, ZERO OUT ALL EDGES BETWEEN NON-CLUSTER NODES
      for cluster in superclusters:
         ADJTMP = []
         for i in range(len(cluster)):
            ADJTMP.append([])
            for j in range(len(cluster)):
               ADJTMP[i].append(self.ADJORIG[self.bacteriaORIG.index(cluster[i])][self.bacteriaORIG.index(cluster[j])])
               
         #print("CLUSTER: "+str(cluster))
         #print("ADJTMP BEFORE:"+str(ADJTMP))
         for i in range(len(ADJTMP)):
            for j in range(len(ADJTMP[0])):
               if (ADJTMP[i][j] < 0):
                  #print ("EDGE SET=0 FROM "+cluster[i]+" TO "+cluster[j])
                  ADJTMP[i][j] = 0  # Remove negative edges

         #print("ADJTMP AFTER: "+str(ADJTMP))
         #x = input()
         subclusters = self.doCluster(ADJTMP, cluster)
         if (len(subclusters) >= 1):
            print("SUPERCLUSTER: "+str(cluster))
            rivalclub = []
            for subcluster in subclusters:
               rivalclub.append(subcluster)
            self.rivalclubs.append(rivalclub)
            print("RIVAL CLUBS: "+str(rivalclub))
            

   def output(self, filename):
      filestuff = open(filename, 'w')
      for rivalclub in self.rivalclubs:
         filestuff.write("\"\",\"x\"\n")
         for subcluster in rivalclub:
            for i in range(0, len(subcluster)):
               filestuff.write(subcluster[i])
               if (i != len(subcluster)-1):
                  filestuff.write(',')
               else:
                  filestuff.write('\n')
