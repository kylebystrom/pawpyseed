import re
import os

directory = 'pawpyseed/core/'

def write_pxd(pxdname, files):
	full_file = ''
	#full_file += '# cython : profile=True\n'
	full_file += '# cython : language_level=3\n'
	full_file += 'from libc.stdio cimport FILE\n'
	for fname in files:
		f = open(directory + fname + '.h', 'r')
		code = f.read()
		code = re.sub('/\*([^*]|[\r\n]|(\*+([^*/]|[\r\n])))*\*+/', '', code)
		code = re.sub('; +', ';', code)
		code = re.sub('(//|#)[^\n]+\n', '', code)
		code = code.replace('{', ':;').replace('typedef', 'ctypedef')
		f.close()
		code_lines = code.split(';')
		code_lines = [c.split('*/')[-1] + '\n' for c in code_lines]
		i = 0

		full_file += '\n\ncdef extern from "%s.h":\n\n\t' % fname

		while i < len(code_lines):
			inc = True
			while code_lines[i].startswith('\n'):
				code_lines[i] = code_lines[i][1:]
			if 'ctypedef' in code_lines[i]:
				linenum = i
			"""
			if '//' in code_lines[i]:
				inc = False
				code_lines[i] = code_lines[i].split('\n', 1)
				if len(code_lines[i]) == 1:
					code_lines.pop(i)
				else:
					code_lines[i] = code_lines[i][1]
			"""
			if '}' in code_lines[i]:
				inc = False
				mystr = code_lines[i].split('}')[-1].replace('\n', '')
				code_lines[linenum] = 'ctypedef struct ' + mystr + ':\n'
				code_lines.pop(i)
			if inc:
				if '(' in code_lines[i]:
					code_lines[i] = 'cdef ' + code_lines[i]
				i += 1

		for line in code_lines:
			full_file += line.replace('\n', '\n\t').replace('(void)', '()')
		
	f = open(os.path.join(directory, pxdname), 'w')
	f.write(full_file.replace('\t', '    '))
	f.close()

write_pxd('pawpyc_extern.pxd',
	['utils', 'projector', 'pseudoprojector', 'reader', 'density', 'sbt', 'linalg', 'radial'])
write_pxd('tests/testc_extern.pxd', ['tests/tests', 'utils'])
