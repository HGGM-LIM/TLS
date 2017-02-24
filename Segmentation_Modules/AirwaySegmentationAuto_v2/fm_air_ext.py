#!/usr/bin/python
import si,os,cPickle,numpy,sys,itk,collections
from si.data.heap import PriorityQueue
from si.data.stack import Slice3
from si.persistence import *
from elixir import session

#from guppy import hpy

from copy import copy

offsets = numpy.array([ [-1,0,0],[1,0,0],[0,0,-1],[0,0,1],[0,-1,0],[0,1,0]])
offsets21 = numpy.array([ [-1,-1,-1],[-1,-1,0],[-1,-1,1],[-1,0,-1],[-1,0,0],[-1,0,1],[-1,1,-1],[-1,1,0],
   [-1,1,1],[0,-1,-1],[0,-1,0],[0,-1,1],[0,0,-1], [0,0,1],[0,1,-1],[0,1,0],[0,1,1],[1,-1,-1],[1,-1,0],
   [1,-1,1],[1,0,-1],[1,0,0],[1,0,1],[1,1,-1],[1,1,0],[1,1,1] ])

currSeg = None
neighs_cache = {}
neighs21_cache = {}


class functiondict(collections.defaultdict):
    def __missing__(self,key):
        if self.default_factory is None: raise KeyError((key,))
        if not (key in self.trials):
            return numpy.inf
        
        value = self.default_factory(key)
        return value


def getNeighs(p):
    """Return 6 connected neighbours""" 
    
    if neighs_cache.has_key(p):
        return neighs_cache[p]
    
    
    maxX = larr.shape[0] - 1
    maxY = larr.shape[1] - 1
    maxZ = larr.shape[2] - 1
    
    pa = numpy.asarray([p])
    neighs = pa.repeat([len(offsets)],axis=0) + offsets
    
    #print "Neighbours",neighs
    
    neighs.ravel()[neighs.ravel() < 0] = 1
    
    neighs[:,0][neighs[:,0] > maxX] = maxX
    neighs[:,1][neighs[:,1] > maxY] = maxY
    neighs[:,2][neighs[:,2] > maxZ] = maxZ
    
  
    
    assert len(neighs[numpy.where(neighs[:,0]<0)]) == 0,"No negative index!"
    assert len(neighs[numpy.where(neighs[:,1]<0)]) == 0,"No negative index!"
    assert len(neighs[numpy.where(neighs[:,2]<0)]) == 0,"No negative index!"
    
    neighs_cache[p] = neighs
    
    return neighs

def getNeighs21(p):
    """Return 21 connected neighbours""" 
    
    
    #lp = tuple(p)
    
    if neighs21_cache.has_key(p):
        return neighs21_cache[p]
    
    maxX = larr.shape[0] - 1
    maxY = larr.shape[1] - 1
    maxZ = larr.shape[2] - 1
    
    pa = numpy.asarray([p])
    neighs = pa.repeat([len(offsets21)],axis=0) + offsets21
    
    neighs.ravel()[neighs.ravel() < 0] = 1
    
    neighs[:,0][neighs[:,0] > maxX] = maxX
    neighs[:,1][neighs[:,1] > maxY] = maxY
    neighs[:,2][neighs[:,2] > maxZ] = maxZ
    
    
    assert len(neighs[numpy.where(neighs[:,0]<0)]) == 0,"No negative index!"
    assert len(neighs[numpy.where(neighs[:,1]<0)]) == 0,"No negative index!"
    assert len(neighs[numpy.where(neighs[:,2]<0)]) == 0,"No negative index!"
    
    neighs21_cache[p] = neighs
    
    return neighs

def calc_viscosity_step(p):
    """ If the point is below a threshold,than return 1, else return infinite """
     
    pixel_val = larr[tuple(p)]
    #print pixel_val
    if (pixel_val < high_threshold) & (pixel_val > low_threshold):
        
        return 1.0
    
    else:
        return numpy.inf

def sol(coeff,v):
    
    #print "Solving potential equation with coeff=",coeff,"p=",v
    
    if numpy.isinf(v):
        return numpy.inf
    
    while (numpy.inf in coeff):
        coeff.remove(numpy.inf)
    
    if len(coeff) == 1:
        # u = Ua + P
        return coeff[0] + v
    
    a = len(coeff)    
    b = -2*sum(coeff)
    c = sum((x**2 for x in coeff)) - v**2
    
    #print "D = b**2 - 4ac with (a,b,c)=",a,b,c
    
    D = b**2 - 4*a*c
    
    if D < 0:
        return None #No real solutions
    
    d = D**(0.5)
        
    #Return the largest solution
    u = max( ( -b - d ), ( -b + d ) ) / (2.0*a)

    return u 

def calc_potential(point,T):
            
    a1,a2,b1,b2,c1,c2 = [tuple(np) for np in getNeighs(point)]
    
    Ua = min(T[a1],T[a2])       
    Ub = min(T[b1],T[b2])       
    Uc = min(T[c1],T[c2])       

    v = calc_viscosity_step(point)

    coeff = [Ua,Ub,Uc]
    coeff.sort()
 
    #print "(Coeff,visc)=",coeff,v
    while(coeff):
        u = sol(coeff,v)

        if u is None:
            del coeff[-1] # Remove last coefficent
            continue

        if u >= coeff[-1]:
            break
        else:
            del coeff[-1] #Remove last coefficent
            continue
    
    return u

def dump_image(segs,storeSegment,spacing,study):
    
    #print "Dumping image",len(segs),storeSegment,spacing
    numPixels = 0 #Count 1 pixels to easy the calculatio of volume later on...
    
    with_segs = numpy.zeros(larr.shape)
    allones = numpy.zeros(larr.shape)
    
    #seg = None
    
    segs_dict = {} 
    if storeSegment:
        for seg in segs:
            airseg = AirwaySegment()
            airseg.dimension = len(seg.alive)*spacing[0]*spacing[1]*spacing[2]/1000
            airseg.hid = seg.id
            airseg.stillToGrow = len(seg.trials)
            airseg.stoppedReason =  seg.stoppedBecause
            airseg.study = study
            
            segs_dict[seg.id] = airseg
    
    for seg in segs:
        if storeSegment:
            currs = segs_dict[seg.id]
            if not (seg.parent is None):
                if not segs_dict.has_key(seg.parent.id):
                    print "Segment",seg.id,'requires',seg.parent.id,'as a parent but it has probably been removed due to a leakage...'
                else:
                    currs.parent = segs_dict[seg.parent.id]
                
        vox = numpy.asarray([tuple(tp) for tp in seg.alive]).astype('i')
        with_segs[vox[:,0],vox[:,1],vox[:,2]] = seg.id
        allones[vox[:,0],vox[:,1],vox[:,2]] = 1
        
 
    with_segs_img = with_segs.astype(numpy.uint16)
    allones_img = allones.astype(numpy.uint16)
    
    with_segs_img = with_segs_img.transpose()
    allones_img = allones_img.transpose()
    
    
    numPixels = allones.sum()
    
    return si.image.fromNumpyArray(with_segs_img),si.image.fromNumpyArray(allones_img),numPixels

class Wavefront(object):
    def __init__(self,pointset):
        self.pointset = pointset
    
    @property
    def radius(self):
        return numpy.sqrt(len(self.pointset)/numpy.pi)
            
    @property
    def centroid(self):
        
        parr = numpy.asarray([p for p in self.pointset]).astype('i')
        
        num = len(parr)
        x = parr[:,0].sum()
        y = parr[:,1].sum()
        z = parr[:,2].sum()
               
        return x / num,y / num, z / num
    
    def _walkConnected(self,point):
   
        visited = set()
  
        nexts = [point]
        finito = False
        
        
        while (not finito):
            
            if len(nexts) <= 0:
                finito = True
                break
            
            p = nexts.pop()
            #print "\tpopped element",p
            
            #print "\tFinding neighborhoods..."
            neighs = getNeighs21(p)
            
            visited.add(p)
        
            #Prendi solo quelli che fanno parte del wavefront
            for n in neighs:
                tn = tuple(n) 
                if (tn in self.pointset) & (not tn in visited):
                    #print "One of the neighs is in the wavefront:",n
                    nexts.append(tn)
            
        return visited
    
    def splitConnectedRegions(self,ignore_smaller=1):
        regions = []
        copy_set = self.pointset.copy()
        
        #print "\tSplitting connected regions:"
        
        while (len(copy_set)>0):
            #print "Remain",len(copy_set),"points to visit"
            str_p = copy_set.pop()
            
            #print "Starting from point",str_p
            
            visited = self._walkConnected(str_p)
            regions.append(visited)
            copy_set.difference_update(visited)
    
        if (ignore_smaller < 1):
            big_regions = []
            small_regions = []
            for r in regions:
                np = len(r)
                perc = (np*1.0 / len(self.pointset))
                #print "\tRegion contains ",np,"points (",perc*100,"%)",Wavefront(r).centroid
                
                if perc > ignore_smaller:
                    big_regions.append(r)
                else:
                    small_regions.append(r)
            
            print "\tReturning #",len(big_regions),"big and #",len(small_regions),"small regions"
            return big_regions,small_regions 
        else:
            return regions,None
            
    
    def is_connected(self):
        #print "Checking connectness"
        p = self.pointset.pop()
        self.pointset.add(p)
        visited = self._walkConnected(p)
        
        return len(visited) == len(self.pointset)

class Segment(object):
    
    def __init__(self,id,seed=None,parent=None,active_wf=None):
        self.leaking = False
        self.id = id
        self.parent = parent
        self.stoppedBecause = ''
        
        self.name = 'Segment'+str(self.id)
        
        if not (parent is None):
            self.currTime = parent.currTime
            self.lastTime = parent.lastTime
            
            self.alive = set()
            self.trials = copy(active_wf.pointset)
            
            #for p in self.trials:
            #    for pn in (tuple(nt) for nt in getNeighs(p)):
            #        if pn in parent.alive:
            #            self.alive.add(pn)
            
            self.T=parent.T        
            self.TTrials=PriorityQueue()
            
            for t in self.trials:
                self.TTrials[t] = self.T[t]
            
            self.lastGrown = set()
            
        else:
            self.currTime = self.lastTime = 0.0
            self.alive = set()
            self.trials = set()
            self.T=functiondict(calc_viscosity_step)
                  
            self.TTrials=PriorityQueue()
            
            self.lastGrown = set()
            self.alive.add(seed)
            self.T[seed] = 0.0
            
            self.trials = set([tuple(tp) for tp in getNeighs(seed)])
            self.trials = self.trials.difference(self.alive)
            
            self.T.trials = self.trials
            
            for t in self.trials:
                self.TTrials[t] = calc_viscosity_step(t)

    
    def setTrials(self,points):
        self.trials = points
        
        copyset = self.TTrials.copy()
        
        self.TTrials.clear()
        
        for t in self.trials:
            if t in copyset:
                self.TTrials[t] = copyset[t]
            else:    
                self.TTrials[t] = calc_viscosity_step(t)
        
        self.T.trials = self.trials
                
       
    
    def removeFromTrials(self,p):
        del self.TTrials[p]
        self.trials.remove(p)
    
 
    
    def growFor(self,numsec):
        """ Grow the wavefront for numsec seconds (timesteps) """
        
        self.lastTime = self.currTime
        
        newTime = self.currTime + numsec
        
        #print "Growing from",self.currTime,"to",newTime
        
        self.lastGrown.clear()
        
        while (self.currTime < newTime):
            
            if len(self.trials) <= 0:
                return True
            
            try:
                minP = self.TTrials.smallest()
            except IndexError:
                print "BUG!: IndexError in priority queue"
                print str(self.TTrials)
                
                            
            self.T[minP] = self.TTrials[minP] 
            
            self.alive.add(minP)       #Froze U(p) by moving p to the alive points
            self.removeFromTrials(minP)
            self.lastGrown.add(minP)
            
            neighs = (tuple(p) for p in getNeighs(minP))
                
            for n in neighs:
                if (n not in self.alive) & (n not in allSegmentedPoints):
        
                    self.trials.add(n)
                   
                    un = calc_potential(n,self.T)
                    #print "U(",n,")=",un
                    self.T[n] = un
                    
                    self.TTrials[n] = un
                    
                    if numpy.isinf(un):
                        self.alive.add(n)
                        self.removeFromTrials(n)
                    else:
                        if un > self.currTime:
                            self.currTime = un    
                    
            
            assert len(self.trials) == len(self.TTrials),"Trials points mismatch:"+str(self.trials)+"\n"+str(self.TTrials)

                
            
    def getWavefront(self):
        #print "Returning wavefront fron time",self.lastTime,"to",self.currTime
                
        wf = Wavefront(self.trials)
        
        
        return wf

def main():
    
    
    #max_gr = 0
    global currSeg,allSegmentedPoints
    
    
        
    print "Starting fast marching level set"
    
    
    
    
    finishedSegments = []
    toShow = []
    
    
    
    
    seg_num = 1
    i = 0
    
    
    seg = Segment(seg_num,initial_seed)
    seg_num += 1                #Inc after creation of a segment
    
    
    segmentQueue = [seg]
    
    while(len(segmentQueue)>0):
        
        #This is quite important: we are using a lifo queue so that early discovered segments (usually the most importants)
        # are processed before the newest ones (often with problems of leaking)
        currSeg = segmentQueue.pop(0)
        
        
        alive_len = len(currSeg.alive)
        trials_len = len(currSeg.trials)
        
        print "START SEGMENT",currSeg.id,'(alive:',alive_len,'trials:',trials_len,')'
        segmentGrown = False
    
        owf_len = -1
        if alive_len > 0:
            aarr = numpy.asarray([tuple(t) for t in currSeg.alive])
            print "alive points start at slice",aarr[:,2].min(),"and stop at slice",aarr[:,2].max()
        
        if trials_len > 0:
            tarr = numpy.asarray([tuple(t) for t in currSeg.trials])
            
            print "trial points start at slice",tarr[:,2].min(),"and stop at slice",tarr[:,2].max()
    
        while(not segmentGrown):
            i += 1
            
            print "> Iteration",i,'Segment',currSeg.id,"trial#",len(currSeg.trials),"alive#",len(currSeg.alive),"segments queue len:",len(segmentQueue)
            
            if i == stopAt:
                segmentGrown = True
                segmentQueue = []
                currSeg.stoppedBecause = 'exceeded user specified max iterations'
                continue
            
            if len(currSeg.trials) == 0:
                segmentGrown = True
                currSeg.stoppedBecause = 'no more points to grow'
                continue
            
            prematurelyEnded = currSeg.growFor(timeStep)
            
            if prematurelyEnded:
                print "Segment couldn't grow for the entire time step" \
                    " because trial points terminated earlier"
                
                currSeg.stoppedBecause = 'no more points to grow (prematurely ended)'
                
                continue
            
                        
             
            wf = currSeg.getWavefront()
            
            wf_len = 1.0*len(wf.pointset)
            
            if owf_len == -1:         #First iteration
                owf_len = wf_len
                                
            wf_grown_perc = (wf_len - owf_len) / owf_len * 100
            
            #print "Wavefront contains",wf_len,"points and grew of",wf_grown_perc,"%"
            
            if (currSeg.currTime > 5.0) & (wf_grown_perc > 150):
                print '*'*5,"Leakage detected:removing offending segment..."
                
                currSeg.leaking = True
                segmentGrown = True
                currSeg.stoppedBecause = 'Leaking'
        
                continue
            
            #Do not check for parent-child growing limit if we are processing early segments (<5)
            if (currSeg.id > 5) & (not (currSeg.parent is None)):
                #grw = len(currSeg.alive) *1.0 / len(currSeg.parent.alive)
                #print grw
                #if grw > max_gr:
                #    max_gr = grw
                if len(currSeg.alive) > 5*len(currSeg.parent.alive):
                    print "Segment is growing more than parent...stopping it"
                    segmentGrown = True
                    currSeg.stoppedBecause = 'Growing more than parent'
                    continue
            
            connected = wf.is_connected()
            
            if not connected:
                
                wavefronts,small_wvs = wf.splitConnectedRegions(ignore_smaller=0.2)
                    
                if (not (small_wvs is None)) & len(small_wvs) > 0:
                    
                    #print "\tSmaller branch will be directly added to the grown segment and removed from the wavefront..."
                 
                    for small in small_wvs:
                        currSeg.alive = currSeg.alive.union(small)
                        currSeg.setTrials(currSeg.trials.difference(small))
            
                if len(wavefronts) == 1:
                    #print "\tEven if the WF is not connected, only one front is to be propagated,"\
                    #     "so I won't start a new segment"
                    pass
                     
                elif (len(wavefronts) > 1):
                    print "Wavefront is NOT connected"
                 
                    print "\tCreating ",len(wavefronts),"new segments..."
                    
                    for w in wavefronts:
                        new_wf = Wavefront(w)
                                                
                        seg = Segment(seg_num,seed=None,parent=currSeg,active_wf=new_wf)
                        seg_num += 1                #Inc after creation of a segment
                        
                        assert len(seg.trials) > 0
                        segmentQueue.append(seg)
                        
                        currSeg.setTrials(currSeg.trials.difference(new_wf.pointset))
                        currSeg.alive = currSeg.alive.union(new_wf.pointset)
                        
                        assert len(seg.trials) > 0
                    
            owf_len = wf_len
                        
            
        finishedSegments.append(currSeg)
        
        allSegmentedPoints = allSegmentedPoints.union(currSeg.alive)
        
    print "Creating segmented image..."
    
    print "Saving",len(finishedSegments),"segments"
    
    for s in finishedSegments:
        if s.leaking:
            finishedSegments.remove(s)
            print "Removed segment",s.id,"because is leaking!"
            continue
        
        if len(s.alive) == 0:
            finishedSegments.remove(s)
            print "Removed segment",s.id,"because has 0 alive points"
            continue
            
    for s in finishedSegments:
        print "Segment",s.id,"alive#",len(s.alive),"trials",len(s.trials)
        
        aarr = numpy.asarray([tuple(t) for t in s.alive])
        print "alive points start at slice",aarr[:,2].min(),"and stop at slice",aarr[:,2].max()
        
        if len(s.trials) > 0:
            tarr = numpy.asarray([tuple(t) for t in s.trials])        
            print "trial points start at slice",tarr[:,2].min(),"and stop at slice",tarr[:,2].max()
    
    return finishedSegments
    

def get_trachea(filename,startFromSlice=0):
    lungs = si.image.asUS(si.image.load3D(filename))
    lungs_arr= si.image.asNumpyArray(lungs)
       
    slice = Slice3(lungs_arr[startFromSlice])
    trachea_point_list = slice.getPointsValued(1) 
    
    npoints = len(trachea_point_list)
    tarr = numpy.asarray([tuple(p) for p in trachea_point_list])
    
    
    return (tarr[:,1].sum()/npoints, tarr[:,0].sum()/npoints, startFromSlice)

    
if __name__ == '__main__':
    if len(sys.argv) < 0:
        print "Usage:",sys.argv[0],"label_file input_file","[stopAtIterationNum=0] [dbfile study_id] "
    else:
        
        neighs21_cache_name = os.path.join(si.properties['cachedir'],'neighs21.dat')
        neighs_cache_name = os.path.join(si.properties['cachedir'],'neighs.dat')
        
        #try:
        neighs21_cache = cPickle.load(file(neighs21_cache_name))
        #except:
        #    pass
        
        #try:
        neighs_cache = cPickle.load(file(neighs_cache_name))
        #except:
        #    pass

        #hp = hpy()
        
        #hp.setrelheap()
        
        
        
        label_file = sys.argv[1]
        input_file = sys.argv[2]
        #initial_seed = tuple((int(c) for c in sys.argv[3].split(',')))
        
        initial_seed = get_trachea(label_file)
        
        #print initial_seed,initial_seed2
        
    
        """ Global variables """
        low_threshold = -2000
        high_threshold = -800
        
        timeStep = 0.5
        
        #wave_front = 48
        stopAt = 0
        try:
            stopAt = int(sys.argv[3])
        except:
            pass
         
        storeSegments = False
        study_id = None
        try:
            dbfile = sys.argv[4]
            study_id = int(sys.argv[5])
            setDbToFile(dbfile)
            setup_all()
            storeSegments = True
        except:
            print "Database storing disabled!"
        
        inputImage = si.image.load3D(input_file)
        spacing = tuple(inputImage.GetSpacing())
        print "Spacing",spacing
        
        #lungs = si.image.asUS(inputImage)
        lungs = inputImage
        larr = si.image.asNumpyArray(lungs) 
        larr = larr.transpose().copy()
    

        allSegmentedPoints = set()
    
        finishedSegments = main()
        
        if (storeSegments):
            study = Study.query.filter_by(id=study_id).one()
        else:
            study = None
        
        segments,allones,numPixels = dump_image(finishedSegments,storeSegments,spacing,study)
        
        if storeSegments:
            session.commit()
        
        itk.write(segments,'air.tiff')
        
        mf = itk.MaskNegatedImageFilter.ISS3IUS3ISS3.New()
        mf.SetInput1(lungs)
        mf.SetInput2(allones)
        
        itk.write(mf.GetOutput(),'lungs-masked-notra.dcm')
        
        
        #f = file('neighs.dat','w')
        #pickle.dump(neighs_cache, f,protocol=2)
        #f.close()
        
        print "Extracted airways volume (ml) is:",spacing[0]*spacing[1]*spacing[2]*numPixels/1000
    
    
