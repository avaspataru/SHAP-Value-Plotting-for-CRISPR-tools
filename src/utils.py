import sys

#returns the object corresponding to the requested tool
def getToolObject(toolName):
    if toolName in ['tuscan-regression','tuscan-classification']:
        sys.path.insert(0,'./tuscan-model')
        from tuscandata import TuscanData
        tool = TuscanData()
        tool.setRegressionFlag(toolName=='tuscan-regression') #if regression
    elif toolName == 'sgrnascorer2':
        sys.path.insert(0, './sgRNAScorer2-model')
        from sgRNAScorer2data import SgRNAScorer2Data
        tool = SgRNAScorer2Data()
    elif toolName == 'wu-crispr':
        sys.path.insert(0,'./wu-crispr-model')
        from wucrisprdata import WuCrisprData
        tool = WuCrisprData()
    elif toolName == 'ssc':
        sys.path.insert(0, './ssc-model')
        from sscdata import SSCData
        tool = SSCData()
    elif toolName in ['chop-chop-xu', 'chop-chop-doench', 'chop-chop-moreno']:
        #ASSUME Xu scoring method = TO DO CHANGE
        sys.path.insert(0, './chop-chop-model')
        from chopchopdata import ChopChopData
        tool = ChopChopData()
        tool.setScoring(toolName)
    else:
        print("Tool Name not valid")
        quit()
    return tool
