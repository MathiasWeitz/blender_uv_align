# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
	'name': "Aligning UV-coords",
	'author': "Mathias Weitz",
	'version': (1, 0, 4),
	'blender': (2, 7, 2),
	'api': 52859,
	'location': "IMAGE_EDITOR > UI",
	'description': "various tricks on UV",
	'category': 'UV'}
	
import bpy 
from bpy.props import *
import bmesh
import math 
import mathutils 
import itertools
from math import pi, cos, sin, sqrt
from mathutils import Vector, Matrix

def debug_del():
	for w in bpy.context.scene.objects:
		if w.name[0:9] == "textdebug":
			bpy.context.scene.objects.unlink(w)
			bpy.data.objects.remove(w)
	
def debug_show(origin,text):
	ft = bpy.data.curves.new('mytext','FONT')
	ft.body = text
	ft.size = 0.1
	ft.shear = 0
	ft.align = 'CENTER'
	ft.resolution_u = 2
	gr = bpy.data.objects.new("textdebug",ft)
	gr.location = Vector(origin)
	gr.show_x_ray = True
	target = bpy.data.objects["Camera"]
	cns = gr.constraints.new("TRACK_TO")
	cns.name = "MTrack"
	cns.target = target
	cns.track_axis = 'TRACK_Z'
	cns.up_axis = 'UP_Y'
	cns.owner_space = 'LOCAL'
	cns.target_space = 'WORLD'
	bpy.context.scene.objects.link(gr)
	
def debug_show_verts(me):
	for v in me.vertices:
		debug_show(v.co,'p' + str(v.index))

def debug_show_polygons(me):
	for v in me.polygons:
		debug_show(v.center,'f' + str(v.index))

class matrix:
	'''my own matrix calculation,
	especially needed Cholesky decomposition and support for wxMaxima'''
	def __init__(self):
		self.field = [[0]]

	def set(self,y,x,value):
		while len(self.field) <= y:
			self.field.append([0])
		while len(self.field[y]) <= x:
			self.field[y].append(0)
		self.field[y][x] = value

	def get(self,y,x):
		erg = 0
		if y < len(self.field):
			if x < len(self.field[y]):
				erg = self.field[y][x]
		return erg

	def dimX(self):
		erg = 0
		for i in range(len(self.field)):
			if erg < len(self.field[i]):
				erg = len(self.field[i])
		return erg

	def copy(self):
		erg = matrix()
		for i in range(len(self.field)):
			for j in range(len(self.field[i])):
				erg.set(i,j,self.field[i][j])
		return erg

	def dimY(self):
		return len(self.field)

	def t_copy(self):
		# transponierte Matrix
		erg = matrix()
		for i in range(len(self.field)):
			for j in range(len(self.field[i])):
				erg.set(j,i,self.field[i][j])
		return erg

	def write(self):
		for v1 in self.field:
			e = ''
			for v2 in v1:
				e += "%7.3f" % v2
			print (e)
			
	def writeWMaxima(self, name='d'):
		follow1 = False
		erg = name + ':matrix('
		for v1 in self.field:
			if follow1:
				erg += ','
			follow1 = True
			erg += '['
			follow2 = False
			dimX = self.dimX()
			for v2i in range(dimX):
				if follow2:
					erg += ','
				follow2 = True
				v2 = 0.0
				if v2i < len(v1):
					v2 = v1[v2i]
				if abs(v2) < 1e-12:
					erg += '0.0'
				else:
					erg += str(v2)
			erg += ']'
		return erg + ');'

	def mul(self,y,a):
		# multipliziert Zeile y mit a
		if y < len(self.field):
			for i in range(len(self.field[y])):
				self.field[y][i] *= a

	def sub(self,y1,y2):
		# zieht von Zeile y1 die Zeile y2 ab
		if y2 < len(self.field):
			maxl = len(self.field[y2])
			if y1 < len(self.field):
				maxl = max(maxl,len(self.field[y2]))
			for i in range(maxl):
				self.set(y1,i, self.get(y1,i) - self.get(y2,i))

	def change(self,y1,y2):
		# tauscht Zeile y1 mit Zeile y2
		for i in range(self.dimX()):
			h = self.get(y1, i)
			self.set(y1, i, self.get(y2, i))
			self.set(y2, i, h)

	def mul_matrix(self, m):
		# multipliziert zwei Matrizen miteinander
		erg = matrix()
		if self.dimX() == m.dimY():
			for y in range(self.dimY()):
				for x in range(m.dimX()):
					s = 0
					for k in range(self.dimX()):
						s += self.get(y,k) * m.get(k,x)
					erg.set(y,x,s)
		return erg

	def mul_vektor(self, v):
		# multipliziert eine Matriz mit einem Array (= Vektor)
		erg = []
		#print ("mul_vektor", self.dimX(), len(v))
		if self.dimX() == len(v):
			for y in range(self.dimY()):
				s = 0
				for k in range(self.dimX()):
					s += self.get(y,k) * v[k]
				erg.append(s)
		return erg

	def getZeros(self):
		# get indeterminable columns
		erg = []
		if self.dimX() == self.dimY():
			for y in range(self.dimY()):
				isZero = True
				for x in range(self.dimY()):
					#print ("zero test", self.get(y,x), self.get(x,y))
					isZero &= abs(self.get(y,x)) < 1e-4
					isZero &= abs(self.get(x,y)) < 1e-4
				if isZero:
					erg.append(y)
		return erg
		
	def appendCol(self, col, neg = False):
		# Spalte anhaengen
		dimX = self.dimX()
		for i in range(len(col)):
			self.set(i,dimX,col[i]);
		return self

	def cholesky(self, e):
		'''decomposition Hermitian, positive-definite matrix
		this is more accurate and faster than gauss'''
		l = matrix()
		for k in range(len(self.field)):
			pivot = self.get(k,k)
			if 0 < pivot: 
				l.set(k,k,sqrt(pivot))
				for i in range(k+1,len(self.field)):
					l1 = self.get(i,k)/l.get(k,k)
					l.set(i,k,l1)
					for j in range(k+1,i+1):
						self.set(i,j,self.get(i,j) - l.get(i,k) * l.get(j,k))
		#logging.info(l.writeWMaxima('chol_l'))
		c = []
		for i in range(len(e)):
			s = e[i]
			for j in range(i):
				s -= l.get(i,j) * c[j]
			if abs(l.get(i,i)) < 1e-12: 
				c.append(0.0)
			else:
				c.append(s / l.get(i,i))
		#logging.info(c)
		x = [0] * len(c)
		for i in range(len(e)-1,-1,-1):
			s = c[i]
			for k in range(i+1, len(e)):
				s += l.get(k,i) * x[k]
			if 1e-12 < abs(l.get(i,i)): 
				x[i] = -s / l.get(i,i)
		#logging.info(x)
		return x

	def gauss(self, e):
		# linear Equitation Gauss
		while len(e) < len(self.field):
			e.append(0)
		# forward step
		i = 0
		lenField = len(self.field)
		while i < lenField-1:
			#print ("forward",i)
			#self.write()
			#print (e)
			#logging.info (self.writeWMaxima('zf' + str(i)))
			pivotE = self.get(i,i)
			#print "Pivot", pivotE
			if 1e-5 < abs(pivotE):
				# alles wie es sein soll
				for j in range(i+1, len(self.field)):
					#print "gauss, forward", i,j
					pivotJ = self.get(j,i)
					#print "Pivot2", pivotJ
					if 1e-5 < abs(pivotJ):
						faktor = pivotE/pivotJ
						self.mul(j, faktor)
						e[j] *= faktor
						e[j] -= e[i]
						self.sub(j,i)
				i += 1
			else:
				erg = 0
				# look for a row to swap
				for j in range(i+1, len(self.field)):
					if i < len(self.field[j]) and 1e-5 < abs(self.field[j][i]):
						erg = j
				if 0 < erg:
					self.change(i,erg)
					h = e[i]
					e[i] = e[erg]
					e[erg] = h
				else:
					# Ausnahme, 0-Zeile
					i += 1
		# test real rank
		rank = len(self.field)-1
		test_ranking = True
		while 0 <= rank and test_ranking:
			for i in range(len(self.field[rank])):
				if 1e-5 < abs(self.field[rank][i]):
					test_ranking = False
			if test_ranking:
				rank -= 1
		#logging.info ('rank: ' + str(rank))
		#if 0 < rank and rank < len(self.field)-1:
		#	subMatrixStart = len(self.field)-rank-1
		#	subMatrix = matrix()
		#	for i1 in range(rank+1):
		#		for i2 in range(rank+1):
		#			subMatrix.set(i1, i2, self.get(i1, i2+subMatrixStart))
		#	#logging.info (subMatrix.writeWMaxima('subm'))
		# backward-Schritt
		lin = False
		for i in range(len(self.field)-1,0,-1):
			#print ("backward",i)
			#self.write()
			#print (e)
			#logging.info (self.writeWMaxima('zb' + str(i)))
			pivotE = self.get(i,i)
			if 1e-5 < abs(pivotE):
				for j in range(0,i):
					#print ("gauss, 2backward", i,j)
					pivotJ = self.get(j,i)
					if 1e-5 < abs(pivotJ):
						faktor = pivotE/pivotJ
						self.mul(j, faktor)
						e[j] *= faktor
						e[j] -= e[i]
						self.sub(j,i)

				pivotK = self.get(i,i)
				if 1e-7 < abs(pivotK):
					e[i] /= pivotK
				else:
					e[i] = 0
					#lin = True
		pivotK = self.get(0,0)
		if 1e-7 < abs(pivotK):
			e[0] /= pivotK
		#else:
		#   lin = True
		#print ("e", e)
		#if lin:
		#   e = None
		return e
	
class rotDir2zy:
	def __init__(self):
		self.rotx, self.roty, self.rotz = 0.0, 0.0, 0.0
	
	def setStartEnd(self, start, end):
		self.start = start
		self.end = end
		self.setDirection(end-start)
		return self
	
	def setMesh(self, mesh):
		self.mesh = mesh
		return self
	
	def setDirection(self, dir):
		self.dir = dir.normalized()
		return self
		
	def calculate(self):
		b = 100
		while 0 < b:
			b -= 1
			m = matrix()
			t = []
			
			m.set(0, 0 , -sin(self.rotz)*cos(self.roty))
			m.set(0, 1 , -cos(self.rotz)*sin(self.roty))
			t.append(cos(self.rotz)*cos(self.roty) - self.dir.x)
			
			m.set(1, 0 , cos(self.rotz))
			m.set(1, 1 , 0)
			t.append(sin(self.rotz) - self.dir.y)
			
			m.set(2, 0 , -sin(self.rotz)*sin(self.roty))
			m.set(2, 1 , cos(self.rotz)*cos(self.roty))
			t.append(cos(self.rotz)*sin(self.roty) - self.dir.z)
			#m.write()
			
			ma = m.t_copy().mul_matrix(m)
			vb = m.t_copy().mul_vektor(t)
			mt = ma.copy().appendCol(vb)
			ee = ma.cholesky(vb)
			del m
			
			#print (abs(ee[0]) + abs(ee[1]))
			if abs(ee[0]) + abs(ee[1]) < 1e-15:
				b = 0
			self.rotz += ee[0]
			self.roty += ee[1]
			#print (b, self.roty, self.rotz, cos(self.rotz)*cos(self.roty), sin(self.rotz), cos(self.rotz)*sin(self.roty), self.dir)
		return self
	
	def setOrtho(self, rotx):
		self.rotx = rotx
		return self
	
	def getElevation(self, vec):
		'''
		−sin(w[2])*sin(w)−sin(w[1])*cos(w[2])*cos(w)
		cos(w[1])*cos(w)
		cos(w[2])*sin(w)−sin(w[1])*sin(w[2])*cos(w)
		'''
		self.calculate()
		orthox = -sin(self.roty) * sin(self.rotx) - sin(self.rotz) * cos(self.roty) * cos(self.rotx)
		orthoy = cos(self.rotz) * cos(self.rotx)
		orthoz = cos(self.roty) * sin(self.rotx) - sin(self.rotz) * sin(self.roty) * cos(self.rotx)
		
		# normal
		normx = self.dir.y * orthoz - self.dir.z * orthoy
		normy = self.dir.z * orthox - self.dir.x * orthoz
		normz = self.dir.x * orthoy - self.dir.y * orthox
		
		# take the startvector as one point on the plane
		px = vec.x - self.start.x
		py = vec.y - self.start.y
		pz = vec.z - self.start.z
		
		# check is all vector are orthogonal, than their dotproduct must be zero
		#print ("direction", self.dir.x, self.dir.y, self.dir.z)
		#print ("ortho", orthox, orthoy, orthoz)
		#print ("norm", normx, normy, normz)
		#print ("direction*ortho", self.dir.x*orthox + self.dir.y*orthoy + self.dir.z*orthoz)
		#print ("direction*norm", self.dir.x*normx + self.dir.y*normy + self.dir.z*normz)
		#print ("ortho*norm", orthox*normx + orthoy*normy + orthoz*normz)
		
		# thats the target
		m = matrix()
		m.set(0,0,self.dir.x)
		m.set(1,0,self.dir.y)
		m.set(2,0,self.dir.z)
		
		m.set(0,1,orthox)
		m.set(1,1,orthoy)
		m.set(2,1,orthoz)
		
		m.set(0,2,normx)
		m.set(1,2,normy)
		m.set(2,2,normz)
		
		solve = m.gauss([px,py,pz])
		# first value is along start-end
		# second is along the plane
		# third is height
		#print ("solve", solve)
		
		# show lines
		if False:
			me = self.mesh
			bm = bmesh.new()
			bm.from_mesh(me)
			startV = Vector((self.start.x, self.start.y, self.start.z))
			ve1 = bm.verts.new(startV)
			startV = startV + Vector((solve[0]*self.dir.x, solve[0]*self.dir.y, solve[0]*self.dir.z))
			ve2 = bm.verts.new(startV)
			startV = startV + Vector((solve[1]*orthox, solve[1]*orthoy, solve[1]*orthoz))
			ve3 = bm.verts.new(startV)
			startV = startV + Vector((solve[2]*normx, solve[2]*normy, solve[2]*normz))
			ve4 = bm.verts.new(startV)
			bm.edges.new((ve1, ve2))
			bm.edges.new((ve2, ve3))
			bm.edges.new((ve3, ve4))
			bm.to_mesh(me)
			bm.free()
		
		return solve

class thickface(object):
	__slots__= "v", "uv", "no", "area", "edge_keys"
	def __init__(self, face, uv_layer, mesh_verts):
		self.v = [mesh_verts[i] for i in face.vertices]
		self.uv = [uv_layer[i].uv for i in face.loop_indices]

		self.no = face.normal
		self.area = face.area
		self.edge_keys = face.edge_keys
		
class MessageOperator(bpy.types.Operator):
	bl_idname = "error.message"
	bl_label = "Message"
	type = StringProperty()
	message = StringProperty()
 
	def execute(self, context):
		self.report({'INFO'}, self.message)
		print(self.message)
		return {'FINISHED'}
 
	def invoke(self, context, event):
		wm = context.window_manager
		return wm.invoke_popup(self, width=400, height=200)
 
	def draw(self, context):
		self.layout.label("Notice")
		row = self.layout.row(align=True)
		row.alignment = 'EXPAND'
		row.prop(self, "message")
		#row = self.layout.split(0.80)
		#row.label("press ok and leave dialog") 
		#row.operator("error.ok")
		
class OkOperator(bpy.types.Operator):
	bl_idname = "error.ok"
	bl_label = "OK"
	def execute(self, context):
		return {'FINISHED'}
	
def twoId(id0,id1):
	if id0 < id1:
		result = str(id0) + ',' + str(id1)
	else:
		result = str(id1) + ',' + str(id0)
	return result
		
def getSubCollection(collection, containing, containingNot = []):
	# results must contain all Elements of containing and None of containingNot
	# [[a,b,c], [a,b,d], [a,c,d], [b,c,d]] and [a,b] => [[a,b,c], [a,b,d]]
	result = []
	for elem in collection:
		add = True
		for elemContaining in containing:
			add &= elemContaining in elem
		for elemContainingNot in containingNot:
			add &= elemContainingNot not in elem
		if add:
			result.append(elem)
	return result

def getUVVert(face, vertIndex):
	'''gets the correspondending UV-Vert'''
	result = 0
	for i in range(len(face.vertices)):
		if face.vertices[i] == vertIndex:
			result = face.loop_indices[i]
	return result
		
def getUVLoops(me):
	'''
		returns dict:
			{
				verts: [] ordered verts
				uv: [[],[]] ordered uv-verts, they are ordered in the same direction as the collection above
						and there be mostly two of them, along each side of the selection
				map: {} maps the uv-index to the vert-index
			}
	'''
	#debug_del()
	#debug_show_verts(me)
	#debug_show_polygons(me)
	result = {}
	result['uv'] = []
	result['map'] = {}
	# a list of all selected Edges as id with the form "lowerVertId,higherVertId"
	edgeSelected = []
	# edgeToFace allows us to find neighboring Faces fast
	edgeToFace = {}
	# and a dict from one selected vert to the next
	nextSelectedVert = {}
	# find the selected edges and prepares an ordering
	for edge in me.edges:
		v0, v1 = edge.vertices[0], edge.vertices[1]
		edgeId = twoId(v0,v1)
		if edge.select:
			edgeSelected.append(edgeId)
			if v0 not in nextSelectedVert:
				nextSelectedVert[v0] = []
			if v1 not in nextSelectedVert:
				nextSelectedVert[v1] = []
			nextSelectedVert[v0].append(v1)
			nextSelectedVert[v1].append(v0)
		edgeToFace[edgeId] = []
	# print ('selected Edges', edgeSelected)
	# print (nextSelectedVert)
	# makes a dict from edges to the adjazent faces
	for face in me.polygons:
		for i in range(len(face.loop_indices)):
			v0, v1 = face.vertices[i], face.vertices[(i + 1) % len(face.loop_indices)]
			edgeId = twoId(v0,v1)
			if edgeId not in edgeToFace:
				edgeToFace[edgeId] = []
			edgeToFace[edgeId].append(face)
	#print ('edgeToFace', edgeToFace)
	
	
	# sort the verts along the selected edge
	sortedSelected = []
	for i,j in nextSelectedVert.items():
		if len(j) == 1:
			sortedSelected.append(i)
			break
	c = 100000
	while 0 < c:
		c -= 1
		next = nextSelectedVert[sortedSelected[-1]]
		nextVert = next[0]
		if nextVert in sortedSelected:
			if 1 < len(next):
				nextVert = next[1]
				if nextVert in sortedSelected:
					# circle
					break
			else:
				# finished
				break
		sortedSelected.append(nextVert)
		
	if 1 < len(sortedSelected):
		# this all just makes sense when the polygon contains at least two edges
		#print ('sorted Verts', sortedSelected)
		
		visitedFaces = []
		#for wayIndex in range(len(sortedSelected)-2):
		#	startVert,nextVert,stopVert = sortedSelected[wayIndex], sortedSelected[wayIndex+1], sortedSelected[wayIndex+2]
		#	idEdgeStart = twoId(startVert,nextVert)
		#	idEdgeStop = twoId(nextVert, stopVert)
		
		idEdgeStart = twoId(sortedSelected[0], sortedSelected[1])
		
		# faces contains both sides of the starting edge
		faces = edgeToFace[idEdgeStart]
		for face in faces:
			visitedFaces.append(face.index)
		
		for face in faces:
			#print ('****')
			# the UV-indices along one side of the selected edges,
			# these can be significant more than the vertices because one vertice can have many uv-verts
			startVert,nextVert = sortedSelected[0], sortedSelected[1]
			startUVVert, nextUVVert = getUVVert(face, startVert), getUVVert(face, nextVert)
			result['map'][startUVVert] = startVert
			result['map'][nextUVVert] = nextVert
			uvVerts = [startUVVert, nextUVVert]
			actFace = face
			#print ('starting Face', actFace.index)
			
			segmentCounter = 0
			b = True
			while b:
				#print ('\n+ segment:', segmentCounter, ', actFace:', actFace.index, ', start, next', startVert, nextVert, ', selectedUV:', uvVerts)
				b = False
				if segmentCounter + 2 < len(sortedSelected):
					startVert, nextVert, targetVert = sortedSelected[segmentCounter], sortedSelected[segmentCounter+1], sortedSelected[segmentCounter+2]
					# get NextFace and all UV-Verts
					for nextEdge, nextFaces in edgeToFace.items():
						ev = nextEdge.split(',')
						ev0,ev1 = int(ev[0]),int(ev[1])
						# the next Face with one end of the edge is the next vert on the selected path
						if (ev0 == nextVert or ev1 == nextVert) and (nextFaces[0] == actFace or nextFaces[1] == actFace):
							# now eleminate the edge that is the selected path
							idSelEdge = twoId(ev0,ev1)
							oppVert = ev0
							if nextVert == ev0:
								oppVert = ev1
							#print ('?edge =', idSelEdge, ', oppVert =' , oppVert , targetVert)
							if idSelEdge in edgeSelected:
								# do not cross the line
								if oppVert == targetVert:
									# but move along the line
									#print ('next Segment - act Face:', actFace.index, ', lim Verts:', ev0, ev1 )
									b = True
									segmentCounter += 1
									uvVert1, uvVert2 = getUVVert(actFace, nextVert), getUVVert(actFace, targetVert)
									if uvVert1 != uvVerts[-1]:
										uvVerts.append(uvVert1)
										result['map'][uvVert1] = nextVert
									if uvVert2 != uvVerts[-1]:
										uvVerts.append(uvVert2)
										result['map'][uvVert2] = targetVert
							else:
								# edge has on point on next Vert
								oppFace = nextFaces[0]
								if actFace.index == nextFaces[0].index:
									oppFace = nextFaces[1]
								# to compare visitedFaces is necessary not to go the same way back
								if oppFace.index not in visitedFaces:
									# move to next face
									#print ('actFace:', actFace.index, ', oppFace:', oppFace.index, ', visitedFaces:', visitedFaces)
									oppVert = ev0
									if nextVert == ev0:
										oppVert = ev1
									# now we have a winner
									b = True
									actFace = oppFace
									visitedFaces.append(actFace.index)
									nextUVVert = getUVVert(actFace, nextVert)
									uvVerts.append(nextUVVert)
									result['map'][nextUVVert] = nextVert
									#print ('nextVerts: ', ev0,ev1, ', oppVert:', oppVert, ', nextVert:', nextVert, ', targetVert:', targetVert )
									if oppVert == targetVert:
										# reached the next segment
										segmentCounter += 1
										targetUVVert = getUVVert(actFace, targetVert)
										uvVerts.append(targetUVVert)
										result['map'][targetUVVert] = targetVert
			result['uv'].append(uvVerts)
	result['verts'] = sortedSelected
	#print ('result', result)
	return result

class UVLineZickZack(bpy.types.Operator):
	'''uv alignment for smothing edges'''
	bl_idname = 'uv.lineralign'
	bl_label = 'lineralign'
	bl_options = {'REGISTER', 'UNDO'}
	
	@classmethod
	def poll(cls, context):
		obj = context.active_object
		return (obj and obj.type == 'MESH')
		#return (obj and obj.type == 'MESH' and bpy.context.scene.tool_settings.mesh_select_mode[1])
		
	def execute(self, context):
		#print ('***************')
		amplitude = float(context.scene.zz_amplitude)
		error = 0
		active = bpy.context.active_object
		bpy.ops.object.mode_set(mode='OBJECT')
		me = active.data
		uv_select_sync = bpy.context.scene.tool_settings.use_uv_select_sync
		
		if not me.uv_textures: 
			me.uv_textures.new()
			
		loop = getUVLoops(me)
		loopVerts = loop['verts']
		#print ('************************************************************')
		print ('selected verts', loopVerts)
		#for loopVert in loopVerts:
		#	print (loopVert, me.vertices[loopVert].co)
		coFirst, coLast = me.vertices[loopVerts[0]].co, me.vertices[loopVerts[-1]].co
		
		# find the best plane by suczessiv guessing
		y = rotDir2zy()
		y.setMesh(me)
		y.setStartEnd(coFirst, coLast)
		med = 0
		scope = 3.2 # close to pi
		for j in range(5): # make 5 iterations
			bestValue, bestMed = 10000, med
			for ii in range(-25, 25):
				vii = med + scope * ii / 25
				y.setOrtho(vii)
				error = 0
				for i in range(1, len(loopVerts)-1):
					solve = y.getElevation(me.vertices[loopVerts[i]].co)
					ee = solve[2] / (abs(solve[1]) + 0.1)
					error = error + ee*ee
					#print ("?", solve)
				#print (vii, error)
				if error < bestValue:
					bestValue = error
					bestMed = vii
			#print ('bestValue', bestValue, bestMed)
			med = bestMed
			scope /= 20
		# bestMed is the optial value
		
		# distance on the mesh, is needed to have a measurement of the distance on the UV-Map
		totalDistance = 0.0
		vertToDist = {loopVerts[0]: 0.0}
		for i in range(0, len(loopVerts)-1):
			dist = (me.vertices[loopVerts[i+1]].co - me.vertices[loopVerts[i]].co).length
			#print ('dist', dist)
			totalDistance += dist
			vertToDist[loopVerts[i+1]] = totalDistance
		#print ('distances', vertToDist)
		
		# elevation
		verts2Elevation = {}
		for i in range(1, len(loopVerts)-1):
			solve = y.getElevation(me.vertices[loopVerts[i]].co)
			h = 0
			# @TODO here comes the sign-switch
			if 0 < solve[1]:
				h = solve[2] / (solve[1] + 0.1)
			else:
				h = solve[2] / (solve[1] - 0.1)
			verts2Elevation[loopVerts[i]] = h
		
		for uvsides in loop['uv']:
			#print ('uvsides', uvsides)
			uvStart, uvEnd = me.uv_layers.active.data[uvsides[0]].uv, me.uv_layers.active.data[uvsides[-1]].uv
			uvlength = (uvEnd - uvStart).length
			# the offsets on the uv-map
			dx, dy = (uvEnd - uvStart).x, (uvEnd - uvStart).y
			# the normalized orthogonal
			hx, hy =dy / uvlength, -dx / uvlength
			# print(uvStart, uvEnd, uvlength)
			# print ('dx, dy', dx, dy)
			for uvVertIndex in uvsides:
				uvVert = me.uv_layers.active.data[uvVertIndex]
				vert = loop['map'][uvVertIndex]
				if vert in verts2Elevation:
					# relative distance of the path on the esh
					relLength = vertToDist[vert] / totalDistance
					# height
					h = amplitude * verts2Elevation[vert]
					# print (uvVertIndex, vert, uvVert)
					# print (vert, relLength, vertToDist[vert], h)
					# absolute position on the UV-map
					dxvu, dyuv = uvStart.x + relLength * dx + h * hx, uvStart.y + relLength * dy + h * hy
					#print(dxvu, dyuv)
					uvVert.uv.x, uvVert.uv.y = dxvu, dyuv
			
		#print('uv_data', me.uv_layers.active.data)
		#print('len uv_data', len(me.uv_layers.active.data))
		
		
		bpy.ops.object.mode_set(mode='EDIT')
		return {'FINISHED'}
	
	
class UVTest(bpy.types.Operator):
	'''straight UV-Line'''
	bl_idname = 'uv.linealign'
	bl_label = 'linealign'
	bl_options = {'REGISTER', 'UNDO'}
	
	@classmethod
	def poll(cls, context):
		obj = context.active_object
		return (obj and obj.type == 'MESH')

	def execute(self, context):
		#print ('***************')
		error = 0
		active = bpy.context.active_object
		bpy.ops.object.mode_set(mode='OBJECT')
		me = active.data
		uv_select_sync = bpy.context.scene.tool_settings.use_uv_select_sync
		
		if not me.uv_textures: 
			me.uv_textures.new()
			
		# getting the edges of a mesh-face aligned to the edge of a uv-face
		# the results are in edge with (vertindex_0, vertindex_1, uv_0, uv_1)
		# there could be more than one entry for vertindex_0, vertindex_1
		uv_layer = me.uv_layers.active.data
		#
		# the marked vertices on the UV-Layer
		markedUV = {}
		#
		# [first vert on mesh, second vert on mesh, first UV-vert index, second UV-vert index]
		# for the same edge on the mesh the first two values are the same, the second two are always different
		edges = []
		#
		# the index of the uv-vert back to the vert in the mesh and it's vector
		# every uv-vert can have one vert on the mesh, but every vert on the mesh can have several uv-verts
		# so it's a one to many
		uvVert2vert = {}
		#
		# adjacent verts in the mesh, but only if they are marked
		vert2vert = {}
		#me_vertices = me.vertices 
		#me_edges = me.edges
		#edges = [e for e in me_edges if me.vertices[e.vertices[0]].select and me.vertices[e.vertices[1]].select]
		for i,face in enumerate(me.polygons):
			for loop_index_0 in range(len(face.loop_indices)):
				# the loop_indeces are the indeces of the UV-verts
				# loop_index_0 and loop_index_1 are UV-indeces
				loop_index_1 = (loop_index_0 + 1) % len(face.loop_indices);
				iv0 = face.vertices[loop_index_0]
				iv1 = face.vertices[loop_index_1]
				uvi0 = loop_index_0
				uvi1 = loop_index_1
				if iv1 < iv0:
					iv0 = face.vertices[loop_index_1]
					iv1 = face.vertices[loop_index_0]
					uvi0 = loop_index_1
					uvi1 = loop_index_0
				
				if face.loop_indices[uvi0] not in uvVert2vert:
					uvVert2vert[face.loop_indices[uvi0]] = [iv0, uv_layer[face.loop_indices[uvi0]].uv.copy()]   
				if face.loop_indices[uvi1] not in uvVert2vert:
					uvVert2vert[face.loop_indices[uvi1]] = [iv1, uv_layer[face.loop_indices[uvi1]].uv.copy()]
					
				select0 = False
				select1 = False
				if uv_select_sync:
					select0, select1 = me.vertices[iv0].select, me.vertices[iv1].select
				else:
					select0, select1 = uv_layer[face.loop_indices[uvi0]].select, uv_layer[face.loop_indices[uvi1]].select
					
				if select0:
					markedUV[face.loop_indices[uvi0]] = iv0
				if select1:
					markedUV[face.loop_indices[uvi1]] = iv1
				if select0 and select1:
					edges.append((iv0, iv1, face.loop_indices[uvi0], face.loop_indices[uvi1]))
					if iv0 not in vert2vert:
						vert2vert[iv0] = []
					if iv1 not in vert2vert:
						vert2vert[iv1] = []
					if iv0 not in vert2vert[iv1]:
						vert2vert[iv1].append(iv0)
					if iv1 not in vert2vert[iv0]:
						vert2vert[iv0].append(iv1)
			
		#print (markedUV)
		#print (edges)
		#print ("uvVert2vert", uvVert2vert)
		#print ("len(uvVert2vert)", len(uvVert2vert))
		#print (vert2vert)
		
		# sorting the verts along the edges
		# start
		vertsOrder = []
		for vi, vin in vert2vert.items():
			# if a vert has only one connected point, this edge is the start
			# the connected point (the first in the array) is the second
			if len(vin) == 1:
				vertsOrder.append(vi)
				vertsOrder.append(vin[0])
				break
		# TODO: if there is no start in the loop of the mesh
		# try to find a start in the verts of the UV-Mesh
		
		# sorting the edges and UV-Edges
		if 0 < len(vertsOrder):
			b = True
			maxc = 10000
			while b and 0 < maxc:
				# appending the edges sukzessiv is simple
				b = False
				maxc -= 1
				v = vert2vert[vertsOrder[-1]]
				vn = v[0]
				if vn == vertsOrder[-2]:
					if 1 < len(v):
						vn = v[1]
					else:
						vn = None
				if vn != None:
					vertsOrder.append(vn)
					b = True
			
			# sorting the UV-Edges is more tricky
			# we have to assume that their coordinates on the UV-map are already separated
			uvEdgeOrder = [[],[]]
			dist = 0.0
			if 1 < len(vertsOrder):
				for i in range(len(vertsOrder) - 1):
					dist += (me.vertices[vertsOrder[i]].co - me.vertices[vertsOrder[i+1]].co).length
					# the second uv-line can remain zero in the case of an open edge
					if (len(uvEdgeOrder[0]) == i + 1 and (len(uvEdgeOrder[1]) == i + 1 or len(uvEdgeOrder[1]) == 0)) \
						or i == 0:
						vi0 = vertsOrder[i]
						vi1 = vertsOrder[i+1]
						for e in edges:
							found = False
							if e[0] == vi0 and e[1] == vi1:
								#print ('->', e);
								uv0, uv1 = e[2], e[3]
								found = True
							if e[0] == vi1 and e[1] == vi0:
								#print ('-<', e);
								uv0, uv1 = e[3], e[2]
								found = True
							if found:
								if i == 0:
									if len(uvEdgeOrder[0]) == 0:
										uvEdgeOrder[0] = [uv0, uv1]
									else:
										uvEdgeOrder[1] = [uv0, uv1]
								else:
									d0 = (uv_layer[uv0].uv - uv_layer[uvEdgeOrder[0][i]].uv).length
									d1 = 100.0
									if 0 < len(uvEdgeOrder[1]):
										d1 = (uv_layer[uv0].uv - uv_layer[uvEdgeOrder[1][i]].uv).length
									#print ("add",e, d0, d1, uvEdgeOrder)
									if abs(d0 - d1) < 1e-8:
										if len(uvEdgeOrder[1]) < len(uvEdgeOrder[0]):
											uvEdgeOrder[1].append(uv1)
										else:
											uvEdgeOrder[0].append(uv1)
									else:
										if d0 < d1:
											uvEdgeOrder[0].append(uv1)
										else:
											uvEdgeOrder[1].append(uv1)
								
							#print (i,uv0.uv, uv1.uv)
					else:
						error = 1
			else:
				error = 2    
			#print (uvEdgeOrder)
		else:
			error = 3
			
		#print ('error', error)
		#for uv in markedUV:
		#    print ("uv", uv, markedUV[uv],uv_layer[uv].uv)
		
		# straighten the edges
		if 0 < error:
			if error == 1:
				bpy.ops.error.message('INVOKE_DEFAULT', message = "Something wrong, maybe different Manifolds")   
			else:
				bpy.ops.error.message('INVOKE_DEFAULT', message = "Something wrong, maybe single line couldn't found")   
		else:
			# ein nonManifold bedeutet, das die Linie nicht an der Kante sitzt
			nonManifold = (0 < len(uvEdgeOrder[1]))
			uv = [0,0]
			uv_old = [0,0]
			w0 = uv_layer[uvEdgeOrder[0][-1]].uv - uv_layer[uvEdgeOrder[0][0]].uv
			if nonManifold:
				w1 = uv_layer[uvEdgeOrder[1][-1]].uv - uv_layer[uvEdgeOrder[1][0]].uv
			d = (me.vertices[vertsOrder[0]].co - me.vertices[vertsOrder[1]].co).length
			for i in range(1,len(vertsOrder) - 1):
				#print ("---" , i, uvEdgeOrder[0][i], uvEdgeOrder[1][i])
				ratiod = d / dist
				d += (me.vertices[vertsOrder[i]].co - me.vertices[vertsOrder[i+1]].co).length
				uv_old[0] = uvVert2vert[uvEdgeOrder[0][i]][1]
				uv[0] = uv_layer[uvEdgeOrder[0][0]].uv + w0*ratiod
				c = 1
				if nonManifold:
					uv_old[1] = uvVert2vert[uvEdgeOrder[1][i]][1]
					uv[1] = uv_layer[uvEdgeOrder[1][0]].uv + w1*ratiod
					c = 2
				for j in range(c):
					for uvi in markedUV:
						#print ('+',uvi, uvVert2vert[uvi][0], uvVert2vert[uvEdgeOrder[j][i]][0], uvVert2vert[uvi][1], uv_old[j])
						if uvVert2vert[uvEdgeOrder[j][i]] == uvVert2vert[uvi] and (uvVert2vert[uvi][1] - uv_old[j]).length < 1e-7:
							#print (uvi, uv_layer[uvi].uv)
							uv_layer[uvi].uv = uv[j]
		#print ("len(markedUV)", len(markedUV))
		
		#for edge in me.edges:
		#	print (edge.vertices, edge.use_seam)
		
		bpy.ops.object.mode_set(mode='EDIT')
		return {'FINISHED'}
	
class UVRound(bpy.types.Operator):
	'''put the UV-coordinates of the selected elements to the closest rounded point'''
	bl_idname = 'uv.round'
	bl_label = 'uvround'
	bl_options = {'REGISTER', 'UNDO'}

	def execute(self, context):
		print ("****")
		#debug_show((0,0,0), 'Textc')
		active = bpy.context.active_object
		bpy.ops.object.mode_set(mode='OBJECT')
		interval = int(context.scene.uv_interval)
		me = active.data
		uv_layer = me.uv_layers.active.data
		uv_select_sync = bpy.context.scene.tool_settings.use_uv_select_sync
		for face in me.polygons:
			for i in range(len(face.loop_indices)):
				iv = face.vertices[i]
				v = me.vertices[iv]
				vselect = v.select
				if not uv_select_sync:
					vselect = uv_layer[face.loop_indices[i]].select    
				if vselect:
					#print (face.index, iv, face.loop_indices[i])
					uv = uv_layer[face.loop_indices[i]].uv
					#print (i,iv, face.loop_indices[i], uv.x, uv.y)
					uv.x = round(interval * uv.x) / interval
					uv.y = round(interval * uv.y) / interval
			
		bpy.ops.object.mode_set(mode='EDIT')
		return {'FINISHED'}

class UVTessellate(bpy.types.Operator):
	'''make a mesh based on the selected UV-faces'''
	bl_idname = 'uv.tessellate'
	bl_label = 'uvtessellate'
	bl_options = {'REGISTER', 'UNDO'}

	def execute(self, context):
		print ("****")
		#debug_show((0,0,0), 'Textc')
		active = bpy.context.active_object
		bpy.ops.object.mode_set(mode='OBJECT')
		interval = int(context.scene.uv_interval)
		hideFaces = context.scene.hideFaces
		me = active.data
		uv_layer = me.uv_layers.active.data
		# tri_elem contains a list of all triangles
		# struct like [[vert_0,vert_1,vert_2],[uv_0,uv_1,uv_2]]
		tri_elem = []
		for face in me.polygons:
			if face.select:
				face.hide=hideFaces
				#print ()
				#print (face.index, len(face.vertices), len(face.loop_indices))
				vectors=[]
				for uvi in face.loop_indices:
					vectors.append(Vector((uv_layer[uvi].uv.x,uv_layer[uvi].uv.y,0.0)))
				tri_tess_list = mathutils.geometry.tessellate_polygon([vectors]) 
				for tri in tri_tess_list:
					#print (tri, uv_layer[tri[0]].uv, uv_layer[tri[1]].uv, uv_layer[tri[2]].uv)
					tri_verts = []
					tri_uv = []
					for i in tri:
						tri_verts.append(face.vertices[i])
						tri_uv.append(face.loop_indices[i])
					tri_elem.append([tri_verts,tri_uv])
		#print ("Tris", tri_elem)
		
		k = int(context.scene.uv_tessellate)
		count = 0
		acc = 5e-6
		vertices = {}
		for i1 in range(k+1):
			#self.report({"INFO"}, str(i1))
			for i2 in range(k+1):
				x = i1 / k
				y = i2 / k
				#print ('********** ',x,y)
				p = False
				for tri in tri_elem:
					v0 = uv_layer[tri[1][0]].uv
					v1 = uv_layer[tri[1][1]].uv
					v2 = uv_layer[tri[1][2]].uv
					w1 = v1-v0
					w2 = v2-v0
					w0 = Vector((x,y)) - v0
					det = w1.x*w2.y - w1.y*w2.x
					detx = w0.x*w2.y - w0.y*w2.x
					dety = w1.x*w0.y - w1.y*w0.x
					if det != 0:
						solu_x = detx / det
						solu_y = dety / det
						if -acc <= solu_x and -acc <= solu_y and solu_x+solu_y <= 1.0+acc:
							#bm = mathutils.geometry.intersect_point_tri_2d(Vector((x,y)),v0,v1,v2)
							#print (solu_x,solu_y,' = ',v0.x,v0.y,',',v1.x,v1.y,',',v2.x,v2.y,',')
							ve0 = me.vertices[tri[0][0]].co
							ve1 = me.vertices[tri[0][1]].co
							ve2 = me.vertices[tri[0][2]].co
							we1 = ve1-ve0
							we2 = ve2-ve0
							p = ve0 + solu_x * we1 + solu_y * we2
				if p != False:
					if i1 not in vertices:
						vertices[i1] = {}
					vertices[i1][i2]=p
					count += 1
		
		print ("count",count)
		start = len(me.vertices)
		me.vertices.add(count)
		c = 0
		for i1 in vertices:
			for i2 in vertices[i1]:
				ve = me.vertices[start+c]
				ve.co = vertices[i1][i2]
				vertices[i1][i2] = ve
				c += 1
		
		bm = bmesh.new()
		bm.from_mesh(me)
		bm.verts.ensure_lookup_table()
		bm.faces.ensure_lookup_table()
		for i1 in vertices.keys():
			if 0 < len(vertices[i1]):
				for i2 in vertices[i1].keys():
					verts = []
					if i1 in vertices and i2 in vertices[i1]:
						verts.append(bm.verts[vertices[i1][i2].index])
					if i1 in vertices and i2+1 in vertices[i1]:
						verts.append(bm.verts[vertices[i1][i2+1].index])
					if i1+1 in vertices and i2+1 in vertices[i1+1]:
						verts.append(bm.verts[vertices[i1+1][i2+1].index])
					if i1+1 in vertices and i2 in vertices[i1+1]:
						verts.append(bm.verts[vertices[i1+1][i2].index])
					if 2 < len(verts):
						#print (verts)
						bm.faces.new(verts)
		bm.to_mesh(me)
		bm.free()
		me.update()
		bpy.ops.object.mode_set(mode='EDIT')
		return {'FINISHED'}
		
class SelectShortest(bpy.types.Operator):
	'''select shortest path between two selected verts, obeys hidden faces'''
	bl_idname = 'uv.selectshortest'
	bl_label = 'selectshortest'
	bl_options = {'REGISTER', 'UNDO'}

	def execute(self, context):
		# print ("****")
		#debug_show((0,0,0), 'Textc')
		active = bpy.context.active_object
		bpy.ops.object.mode_set(mode='OBJECT')
		interval = int(context.scene.uv_interval)
		hideFaces = context.scene.hideFaces
		uv_select_sync = bpy.context.scene.tool_settings.use_uv_select_sync
		# print ('uv_select_sync',uv_select_sync)
		me = active.data
		uv_layer = me.uv_layers.active.data
		vert2vert = {}
		vert2uv = {}
		# determine length from one vert to its neighbours
		for i in range(len(me.edges)):
			edge = me.edges[i]
			edge_visible = False
			if uv_select_sync:
				edge_visible = not edge.hide
			else:
				edge_visible = active.data.vertices[edge.vertices[0]].select and active.data.vertices[edge.vertices[1]].select
			if edge_visible:
				#print (edge.vertices[0])
				v0 = active.data.vertices[edge.vertices[0]]
				v1 = active.data.vertices[edge.vertices[1]]
				d = (v0.co - v1.co).length
			
				if v0.index not in vert2vert:
					vert2vert[v0.index] = {}
				if v1.index not in vert2vert[v0.index]:
					vert2vert[v0.index][v1.index] = d           
				if v1.index not in vert2vert:
					vert2vert[v1.index] = {}
				if v0.index not in vert2vert[v1.index]:
					vert2vert[v1.index][v0.index] = d
				
		verts = {}
		start=0
		second=0
		selected = []
		for face in me.polygons:
			for i in range(len(face.loop_indices)):
				iv = face.vertices[i]
				vert = me.vertices[iv]
				vertSelect = vert.select
				if not uv_select_sync:
				   vertSelect = uv_layer[face.loop_indices[i]].select   
				if vertSelect:
					if iv not in selected:
						selected.append(iv)
				if iv not in vert2uv:
					vert2uv[iv] = []
				vert2uv[iv].append(face.loop_indices[i])    
		print ('selected', selected)
		perm = itertools.permutations(selected)
		shortest = {}
		
		if 1 < len(selected) and len(selected) < 9:
			for iv in selected:
				for i in range(len(me.vertices)):
					vert = me.vertices[i]
					if not vert.hide:
						verts[vert.index] = {'d': 10000.0, 'path':[]}
				verts[iv]['d'] = 0.0
				search = [iv]
				c = 1500000
				while 0 < len(search) and 0<c:
					c -= 1
					next = search.pop(0)
					d = verts[next]['d']
					#print ('***', next)
					for nvert in vert2vert[next].keys():
						dp=vert2vert[next][nvert]
						if d + dp < verts[nvert]['d']:
							search.append(nvert)
							verts[nvert]['d'] = d + dp
							verts[nvert]['path'] = [next] + verts[next]['path']
				for iv2 in selected:
					if iv != iv2:
						shortest[str(iv) + '-' + str(iv2)] = verts[iv2]
			bestp = () 
			bestd = 100000.0
			for p in perm:
				d = 0.0
				for i in range(len(p)-1):
					d += shortest[str(p[i]) + '-' + str(p[i+1])]['d']
				#print (p,d)
				if d < bestd:
					bestd = d
					bestp = p[:]
			print ('best', bestp, bestd)
			
			for i in range(len(bestp)-1):
				pp = shortest[str(bestp[i]) + '-' + str(bestp[i+1])]['path']
				for iv in pp:
					if uv_select_sync:
						me.vertices[iv].select = True
					else:
						for iuv in vert2uv[iv]:
							uv_layer[iuv].select = True

		#for i in range(len(me.vertices)):
		#    vert = me.vertices[i]
		#    if not vert.hide:
		#        verts[vert.index] = {'d': 10000.0, 'path':[]}  
		#        if vert.select:
		#            second = start
		#            start = i
		#if 0 < start and 0 < second:
		#    search = [start]
		#    c = 500000
		#    verts[start]['d'] = 0.0
		#    while 0 < len(search) and 0<c:
		#        c -= 1
		#        next = search.pop(0)
		#        d = verts[next]['d']
		#        #print ('***', next)
		#        for nvert in vert2vert[next].keys():
		#            dp=vert2vert[next][nvert]
		#            if d + dp < verts[nvert]['d']:
		#                search.append(nvert)
		#                verts[nvert]['d'] = d + dp
		#                verts[nvert]['path'] = [next] + verts[next]['path']
		#        #print (search)
		#    for v in verts[second]['path']:
		#        #print (v)
		#        me.vertices[v].select = True
			
		#print (vert2vert)
		#print (verts)
		me.update()    
		bpy.ops.object.mode_set(mode='EDIT')
		return {'FINISHED'}          

class VIEW3D_PT_tools_UVTest(bpy.types.Panel):
	bl_space_type = 'IMAGE_EDITOR'
	bl_region_type = 'UI'
	bl_idname = 'uv_even'

	bl_label = "uv align raster"
	bl_context = "objectmode"
	bl_options = {'DEFAULT_CLOSED'}

	def draw(self, context):
		active_obj = context.active_object
		layout = self.layout

		colm = layout.column(align=True)
		row = colm.split(0.25)
		#row.split = 0.15
		w = row.prop(context.scene, "uv_interval")
		#row.split(percentage=0.15)
		#w.alignment = 'RIGHT'
		row.operator("uv.round", text="Snap Selected")
		
		colm = layout.column(align=True)
		col = colm.column(align=True)
		col.operator("uv.selectshortest", text="Select Shortest")

		colm = layout.column(align=True)
		col1 = colm.column(align=True)
		col1.operator("uv.linealign", text="Line Straight")
		colzz2 = colm.column(align=True)
		rowzz2 = colzz2.split(0.4)
		rowzz2.prop(context.scene, "zz_amplitude")
		rowzz2.operator("uv.lineralign", text="Line ZZ")
		
		#col = colm.column(align=True)
		#col.operator("uv.linealign", text="Round")
		
		colm1 = layout.column(align=True)
		colm2 = colm1.column(align=True)
		row2 = colm2.split(0.25)
		w = row2.prop(context.scene, "uv_tessellate")
		row2.operator("uv.tessellate", text="Remesh from UV")
		
		colm3 = colm1.column(align=True)
		row3 = colm3.row(align=True)
		row3.prop(context.scene, "hideFaces")
		

classes = [MessageOperator, OkOperator,
	UVTest,
	UVLineZickZack,
	UVRound,
	UVTessellate,
	SelectShortest,
	VIEW3D_PT_tools_UVTest]   
					
def register():
	#bpy.utils.register_module(__name__)
	for c in classes:
		bpy.utils.register_class(c)
	bpy.types.Scene.uv_tessellate = EnumProperty(
		name="",
		description="mesh size on UV-Editor",
		items=[("256","1","1"),
			   ("128","2","2"),
			   ("64","4","4"),
			   ("32","8","8"),
			   ("16","16","16"),
			   ("8","32","32"),
			  ],
		default='32')
	bpy.types.Scene.uv_interval = EnumProperty(
		name="",
		description="snap factor",
		items=[("256","1","1"),
			   ("128","2","2"),
			   ("64","4","4"),
			   ("32","8","8"),
			   ("16","16","16"),
			   ("8","32","32"),
			  ],
		default="8")
	bpy.types.Scene.zz_amplitude = FloatProperty(
		name="",
		description="amplitude of the zz-line",
		step = 10,
		default = 0.1,
		min = -1.0,
		max = 1.0)
	bpy.types.Scene.hideFaces = BoolProperty(
			name="Hide original Faces",
			description="Hide the selected Faces",
			default=True,
			)
   
def unregister():
	#bpy.utils.unregister_module(__name__)
	for c in classes:
		bpy.utils.unregister_class(c)
   
if __name__ == "__main__":
	register()   