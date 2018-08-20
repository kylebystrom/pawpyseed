import re
import sys

filename = sys.argv[1]

p = re.compile("#[0-9A-Fa-f]{6}")

f = open(filename, 'r')
string = f.read()
f.close()
res = p.findall(string)

mapping = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F',]
mapping = [item for item in zip(mapping, reversed(mapping))]
mapping = {k:v for k,v in mapping}
mapping2 = {'f':'0','e':'1','d':'2','c':'3','b':'4','a':'5'}
for key in mapping2:
	mapping[key] = mapping2[key]
print(mapping)

def map(string):
	mystr = string[1:]
	mystr = [mapping[char] for char in mystr]
	fstr = '#'
	for char in mystr:
		fstr += char
	return fstr

for item in res:
	string = string.replace(item, map(item))
string = string.replace('color: white', 'color: #000000')
string = string.replace('color: black', 'color: #32CD32')

print(string)
f = open(filename, 'w')
f.write(string)
f.close()
"""
p = re.compile("#[0-9a-f]{6}")

f = open(filename, 'r')
string = f.read()
f.close()
res = p.findall(string)

mapping = ['0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f']
mapping = [item for item in zip(mapping, reversed(mapping))]
mapping = {k:v for k,v in mapping}
print(mapping)

def map(string):
	mystr = string[1:]
	mystr = [mapping[char] for char in mystr]
	fstr = '#'
	for char in mystr:
		fstr += char
	return fstr

for item in res:
	string = string.replace(item, map(item))

print(string)
f = open(filename, 'w')
f.write(string)
f.close()
"""
