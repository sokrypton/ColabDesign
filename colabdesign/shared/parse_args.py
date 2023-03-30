import sys, getopt
# class for parsing arguments
class parse_args:
  def __init__(self):
    self.long,self.short = [],[]
    self.info,self.help = [],[]

  def txt(self,help):
    self.help.append(["txt",help])

  def add(self, arg, default, type, help=None):
    self.long.append(arg[0])
    key = arg[0].replace("=","")
    self.info.append({"key":key, "type":type,
                      "value":default, "arg":[f"--{key}"]})
    if len(arg) == 2:
      self.short.append(arg[1])
      s_key = arg[1].replace(":","")
      self.info[-1]["arg"].append(f"-{s_key}")
    if help is not None:
      self.help.append(["opt",[arg,help]])

  def parse(self,argv):
    for opt, arg in getopt.getopt(argv,"".join(self.short),self.long)[0]:
      for x in self.info:
        if opt in x["arg"]:
          if x["type"] is None: x["value"] = (x["value"] == False)
          else: x["value"] = x["type"](arg)

    opts = {x["key"]:x["value"] for x in self.info}
    print(str(opts).replace(" ",""))
    return dict2obj(opts)

  def usage(self, err):
    for type,info in self.help:
      if type == "txt": print(info)
      if type == "opt":
        arg, helps = info
        help = helps[0]
        if len(arg) == 1: print("--%-15s : %s" % (arg[0],help))
        if len(arg) == 2: print("--%-10s -%-3s : %s" % (arg[0],arg[1].replace(":",""),help))
        for help in helps[1:]: print("%19s %s" % ("",help))
    print(f"< {err} >")
    print(" "+"-"*(len(err)+2))
    print("        \   ^__^               ")
    print("         \  (oo)\_______       ")
    print("            (__)\       )\/\   ")
    print("                ||----w |      ")
    print("                ||     ||      ")
    sys.exit()

class dict2obj():
  def __init__(self, dictionary):
    for key in dictionary:
      setattr(self, key, dictionary[key])