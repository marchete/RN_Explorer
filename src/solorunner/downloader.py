import sys
import os
import zlib
import base64
#for downloading level 931+
print("first_level")
width, height = [int(i) for i in input().split()]
valor=str(width)+" "+str(height)+os.linesep
for i in range(height):
	valor+=input()+os.linesep
cmpstr = zlib.compress(valor.encode('utf-8'))
encoded = base64.b64encode(cmpstr).decode("utf-8")
print(encoded, file=sys.stderr, flush=True)
#decoded=base64.b64decode(encoded.encode('utf-8'))
#uncmpstr=zlib.decompress(decoded).decode("utf-8")
#print(uncmpstr, file=sys.stderr, flush=True)