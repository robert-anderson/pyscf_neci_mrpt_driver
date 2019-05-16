

class FileBlock(list):
	def __init__(self):
		list.__init__(self)

	def render(self, fname, **kwargs):
		tmp = '\n'.join([i.format(**kwargs) for i in self])
		if fname is not None:
			with open(fname, 'w') as f: f.write(tmp)
		else: return tmp

	def prefix_all(self, string):
		for i in range(len(self)): self[i] = string+self[i]

	def suffix_all(self, string):
		for i in len(self): self[i] = self[i]+string

	def __lshift__(self, line):
		if isinstance(line, list): self.extend(line)
		else: self.append(line)
		return self

