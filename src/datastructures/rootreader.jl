import UpROOT: UpROOT.uproot
import DataFrames: DataFrame

"""
rootreader(fileName::String, treeName::String)

Read a single ROOT TTree from a single TFile and return
the data as a DataFrame; preserving the names of the variables
and order of the events.
"""
function rootreader(fileName::String, treeName::String)
  tfile = uproot.open(fileName)
  ttree = tfile.get(treeName)
  branches = ttree.keys()
  events = ttree.arrays(branches)
  DataFrame(events)
end
