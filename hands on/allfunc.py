##TODO
"""
grouped_basisの作成は、複雑なのは自作してもらうかto be continued, 金属種ごとくらいなら実装する?金属間化合物や酸化物で重要になるが、今の間はしない。後ほど行う。
サイトの違いはgrouped_basisで対応できる
ひとまずCu-Coに速く映りたい　Cu111から　七種君にもやることを確認。急いで毎日やる
logger infoとdebugを書いていく 最後のチェックがてら。読みにくいところは修正の余地ありなのでそのチェックにも
"""

import ast

def check_type(Object, Type):
    if Object is None:
        raise ValueError("The object is None. Parsing may be failed.") ##もう少し考える
    if Type is list or Type is tuple or Type is dict:
        return ast.literal_eval(Object)
    elif isinstance(Object,Type):
        return Object
    else:
        return Type(Object)


import pprint
import configparser

class Prepare_Config():
	def __init__(self,config_path=None):
		self.config = configparser.ConfigParser()
		if config_path is None:
			self.config.read(os.path.dirname(__file__)+'/config.ini', encoding='utf-8') ##もう少しましな実装にする
		else:
			self.config.read(config_path, encoding='utf-8') ##自作を使いたい人向け。後で引数の説明を書く。

	def __call__(self,section, index, Type):
		parameter = self.config.get(section,index)
		Object = check_type(parameter,Type)
		return Object

	def read(self):
		contents = {}
		for section_name, section_content in dict(self.config).items():
			contents[section_name] = dict(section_content)
		pprint.pprint(contents)

	def save(self,save_path):
		with open(save_path, "w") as f:
			config.write(f)

	def write(self,section,key,value):
		if self.config is None:
			raise NotImplementedError("Overwriting config.ini of SurfQit is prohibited. Please save it in your path, after then write the value in the saved one.")
		self.config[section][key] = str(value)

import os
import ase.io

def storage_path_initialize(config):
	"""
	input:
		config: Prepare_Configクラスのもの。
	description:
		storage_pathを設定し、hostの構造と使用するconfigを保存する(ある意味、その時々の状況のセーブになっている)
	"""
	storage_path = config('preparation','storage_path',str)
	host_path    = config('preparation','host_path',str)
	model_type   = config('preparation','model_type',str)
	if os.path.isdir(storage_path):
		raise FileExistsError("'work' directory exist. Please rename, delete or move it.")
	os.mkdir(storage_path)
	host = ase.io.read(host_path) ##ここができないなら、自分でworkにcifで入れるように警告?
	ase.io.write(storage_path+'host.cif',host,format='cif')
	with open(storage_path+'config.ini', 'w') as file:
		config.config.write(file)

##primから取ってきてQUBO作るなら不要
import spglib
import ase.io
import json

def output_symmetry(config,storage_path=None,host=None):
	"""
	input:
		config: Prepare_Configクラスのもの。
		storage_path: configで自動設定されるが、別用途で使いたいときのため
		host: 同上, 設定の時はase.Atoms class
	output:
		rotations: 回転対称操作
		translations: 並進対称操作
	description:
		host.cifの対称性を確認するとともに対称操作を返す
	"""
	if storage_path==None:
		storage_path   = config('preparation','storage_path',str)
	symmetry_tolerance = config('preparation','symmetry_tolerance',float)
	space_group        = config('preparation','space_group',int)
	model_type         = config('preparation','model_type',str)
	if host==None:
		host      = ase.io.read(storage_path+'host.cif',format='cif')
	model_type           = config('preparation','model_type',str)
	if model_type != "bulk": # slabではz方向の対称性を失わせる
		max_z = np.max(host.get_positions()[:,2])
		for atom in host:
			if atom.position[2] == max_z:
				atom.position[2] += symmetry_tolerance*2 ##いらんかも ################################
	symmetry      = spglib.get_symmetry_dataset((host.get_cell(),np.round(host.get_scaled_positions(), decimals=-1*int(str(symmetry_tolerance).split("e")[1]))%1.0,host.get_atomic_numbers()),symprec=symmetry_tolerance,hall_number=space_group)
	symmetry_data = {}
	symmetry_data['space group symbol'] = str(symmetry['international'])
	symmetry_data['space group number'] = int(symmetry['number'])
	symmetry_data['space group hall']   = str(symmetry['hall'])
	symmetry_data['equivalent atom']    = symmetry['equivalent_atoms'].tolist()
	symmetry_data['origin']             = list(range(len(host)))
	symmetry_data['rotation_operation'] = symmetry['rotations'].tolist()
	symmetry_data['translate_operation']= symmetry['translations'].tolist()
	symmetry_data['equivalent_arrange'] = symmetry_equivalent_arrange(host, symmetry).tolist()
	with open(storage_path+'symmetry.json','w') as f:
		json.dump(symmetry_data, f, indent=2)
	return symmetry['rotations'],symmetry['translations']

import ase.io

def convert_nanosheet_into_slab(config,host=None): ##この関数をoutput_symmetryに組み込むか、組み込むとconfigで設定するかを考える必要がある。が、後で何とでもなるのでとりあえず外に置いておく。また、aseで作った奴でひっくり返すだけで元の位置に重なるものではない限り不要。判定する関数?
	"""
	input:
		config: Prepare_Configクラスのもの。
		host: configで自動設定されるが、別用途で使いたいときのため。設定の時はase.Atoms class
	output:
		None
	description:
		ナノシートモデル(上下とも緩和)をスラブモデル(固定層を設定)に変換。output_symmetryに入れる前に行うことでスラブモデルの場合はz方向に少し原子を動かすことでz方向の対称性をなくす。
	"""
	storage_path   = config('preparation','storage_path',str)
	symmetry_tolerance   = config('preparation','symmetry_tolerance',float)
	if host==None:
		host = ase.io.read(storage_path+'host.cif',format='cif')
	max_z = np.max(host.get_positions()[:,2])
	for atom in host:
		if atom.position[2] == max_z:
			atom.position[2] += symmetry_tolerance*10
	ase.io.write(storage_path+'host.cif',host,format='cif')

import numpy as np
import cython ##このインポートは不要, ただ要件に書く必要がある。多分pipの時に各自の環境に合わせてコンパイルさせないといけない。QUBO作成などもCythonizeした方がいい, https://qiita.com/kenmaro/items/877b47fc293dafc395daによるとsetup.pyでいいらしい。だめならまた考える
##Cythonizeをimportでやるかsetupでやるかは必要。https://hack.nikkei.com/blog/advent20211225/#pyx%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB%E3%81%AE%E3%82%B3%E3%83%B3%E3%83%91%E3%82%A4%E3%83%AB%E6%96%B9%E6%B3%95
##setup.pyの時にbuildする
from cython_process import C_position_to_index 

def symmetry_equivalent_arrange(host, symmetry): # memory-friendly imprementation
	"""
	input:
		host: 対称性を明らかにしたい構造。
		symmetry: spglib.get_symmetry_datasetの結果か、dict:{"rotations": output_symmetry(config)[0],"translations": output_symmetry(config)[1]}でOK
	output:
		symmtetry_equivalent_arrange_index: 並進回転対象運動後の原子位置の移動先
	description:
		host.cifの対称性を確認するとともに対称操作を返す
	"""
	rotate         = symmetry['rotations']    # [operation,xyz,xyz]
	translate      = symmetry['translations'] # [operation,abc]
	fractional     = host.get_scaled_positions()
	operated       = np.dot(rotate,fractional.T).transpose((2,0,1))+translate[np.newaxis,:,:] # [atom,operation,xyz]
	fractional     = np.round(fractional, decimals=10)%1.0 # 丸め誤差で-1e-16とかが出てくると%1.0が1.0になってしまう
	operated       = np.round(operated, decimals=10)%1.0
	operated_index = np.unique(C_position_to_index(fractional,operated.transpose(1,0,2)),axis=1) #[atom, operation]
	return operated_index.T # [operation, (index of ) atoms] ##ここ変えたので、QUBO作成で困るかも。だめなら戻す

import ase.io
import os
#import output_symmetry ##同じファイル内なので今は不要
#import symmetry_equivalent_arrange ##同じファイル内なので今は不要

##ネーミングセンスがない。検討。グラフみたいに思えてしまう
def visualize_symmetry_operation(config):
	"""
	input:
		config: Prepare_Configクラスのもの。
	output:
		None
	description:
		visualizeフォルダにsymmetry_operationフォルダを作り、対称操作の結果を記したcifファイルを保存する
	"""
	storage_path          = config('preparation','storage_path',str)
	roration, translation = output_symmetry(config)
	host = ase.io.read(storage_path+'host.cif',format='cif')
	symmtetry_equivalent_arrange_index = symmetry_equivalent_arrange(host,{"rotations": roration,"translations": translation})
	if not os.path.isdir(storage_path+"visualize/"):
		os.mkdir(storage_path+"visualize/")
	if os.path.isdir(storage_path+"visualize/symmetry_operation/"):
		raise FileExistsError(f"{storage_path}/visualize/symmetry_operation/ directory exist. Please rename, delete or move it.")
	else: 
		os.mkdir(storage_path+"visualize/symmetry_operation/")
	surface = np.where(host.positions[:,2]==np.max(host.positions[:,2]))
	if len(surface)<3:
		surface = np.where(host.positions[:,2]>=np.max(host.positions[:,2])-1)
	visual_preference_surface =np.argsort(np.sum((host.positions[surface][:,:2]-np.array([0,0]))**2,axis=1))
	substitute_indice = [surface[0][(np.where(visual_preference_surface==preference)[0][0])] for preference in [0,1,2]]	
	for model_index,operated_index in enumerate(symmtetry_equivalent_arrange_index): ##ここ変えたので、QUBO作成で困るかも。だめなら戻す
		model = host.copy()
		for atomic_number_difference, substitute_index in enumerate(operated_index[substitute_indice]):
			model[substitute_index].number = host.numbers[substitute_index]+1+atomic_number_difference
		ase.io.write(storage_path+f"visualize/symmetry_operation/symmetry_model{model_index}.cif",model)

from ase.build import surface, add_vacuum

def bulk2slab(conventional_cell,miller,size,vacuum_thickness):
	##primitive atoms object of surface must be primitive along x-y but target structure along z
	##surfaceはprimitiveに対して働かないらしいが、ちゃんと動いているので今は放置。必要があれば(ある意味セルを再構築しているだけなので)自分で書くのもあり
	prim = surface(conventional_cell, miller,size[2],periodic=True)
	prim.set_tags(np.arange(size[2]).repeat(len(conventional_cell))) ## 入力が2*2とかの時にrepeatでいいのか不明, そもそも複数サイトで働かない?⇒basis groupを考える必要。, 酸化物表面も重要。
	prim = wrap_and_sort_by_position(prim)
	size_x,size_y,_ = size
	size = (size_x,size_y,1)
	add_vacuum(prim,vacuum_thickness)
	return prim,size

import numpy as np

def adjust_concentration_each_layer(prim,basis_elements,grouped_basis,ranges):
	##primのtagは層数に対応している。0の中のindexを拡張していくイメージで		
	layer_num = len(set(prim.get_tags()))
	unit_atom_num = len(prim)//layer_num
	if unit_atom_num*layer_num != len(prim):
		raise NotImplementedError("SurfQit does not support surface model with heterogeneous number of atom in each layer")
	basis_elements_each_layer,grouped_basis_each_layer,ranges_each_layer = np.zeros(len(prim)).tolist(),[],[]
	for elements,index,concentrain_range in zip(basis_elements,grouped_basis,ranges): # 同じ組成の層ごとに処理
		for each_index in index: #それぞれの層に元素候補を代入(層ごと)
			basis_elements_each_layer[each_index] = elements
		grouped_basis_each_layer.append(index) # 同じ組成の層ごとに代入
		ranges_each_layer.append(concentrain_range) # 同じ組成の層ごとに代入
		ref_index = index
		"""大本でprimをconventionalの代わりに使うならここを省略できる。allだけ0層目のことを他の層でも繰り返しているのが謎だった
		for layer in range(layer_num)[1:]: # 拡大した場合のことをやっているのか????不明
			print("230",layer,tuple([ind+unit_atom_num*layer for ind in index]),basis_elements_each_layer)
			objective_index = tuple([ind+unit_atom_num*layer for ind in index]) 
			for each_index in objective_index:
				basis_elements_each_layer[each_index] = elements
			grouped_basis_each_layer.append(objective_index)
			ranges_each_layer.append(concentrain_range)"""
	return basis_elements_each_layer,grouped_basis_each_layer,ranges_each_layer

from clease.settings.atoms_manager import AtomsManager

class SurfaceAtomsManager(AtomsManager):
	def __init__(self,atoms):
		super().__init__(atoms)

	def index_by_basis(self,grouped_basis):
		ind_by_symbol = [list(index) for index in grouped_basis]
		return ind_by_symbol

	def basisindex_by_tag(self,tag2grouped_basis):
		ind_by_tag = [[] for tag_number in range(max(list(tag2grouped_basis.values()))+1)]
		for ind,atom in enumerate(self.atoms):
			ind_by_tag[tag2grouped_basis[atom.tag]].append(ind)
		return ind_by_tag

from copy import deepcopy
from clease.settings import ClusterExpansionSettings
from clease.settings.settings import _get_concentration
from clease.settings.template_atoms import TemplateAtoms
from clease.settings.template_filters import ValidConcentrationFilter
from clease.settings.atoms_manager import AtomsManager
from clease.tools import wrap_and_sort_by_position
#import SurfaceAtomsManager
#import ClusterManager_each_sublattice

class ClusterExpansionSettings_nowarning(ClusterExpansionSettings):
	def __init__(self,prim,concentration,size,supercell_factor = 27,db_name = "clease.db",max_cluster_dia=(5.0, 5.0, 5.0),include_background_atoms=False,basis_func_type="polynomial"):
		self._include_background_atoms = include_background_atoms
		self._cluster_mng = None
		self._trans_matrix = None
		self._cluster_list = None
		self._basis_func_type = None
		self._prim_cell = None
		self.concentration = _get_concentration(concentration)

		self.basis_elements = deepcopy(self.concentration.basis_elements)
		self._check_first_elements()

		self.db_name = db_name
		self._set_prim_cell(prim)

		#index_by_symbolを変えるほうが早そうだが、影響やいかに。。。
		prim_mng = SurfaceAtomsManager(self.prim_cell)
		grouped_basis = deepcopy(concentration.grouped_basis)
		self.tag2grouped_basis(self.prim_cell,concentration.grouped_basis)
		prim_ind_by_basis = prim_mng.index_by_basis(grouped_basis) ##modified
		conc_filter = ValidConcentrationFilter(concentration, prim_ind_by_basis)

		self.atoms_mng = SurfaceAtomsManager(None)

		# Construct the template atoms, this is a protected property
		# Access and changes to size and/or supercell factor are redirected to
		# this instance.
		self._template_atoms = TemplateAtoms(
			self.prim_cell,
			supercell_factor=supercell_factor,
			size=size,
			skew_threshold=40,
			filters=[conc_filter],
		)

		self.set_active_template()

		self.max_cluster_dia = max_cluster_dia

		self.basis_func_type = basis_func_type

		if len(self.basis_elements) != self.num_basis:
			raise ValueError("list of elements is needed for each basis")

		# For storing the settings from the CLEASE factories.
		# Just for bookkeeping
		self.kwargs = {}

	def _check_first_elements(self):
		pass

	def tag2grouped_basis(self,prim_cell,grouped_basis): ##ここでindexごちゃりそうで怖い
		self.tag2grouped_basis = {}
		for group,content in enumerate(grouped_basis): ##ここのタグづけがcf_calculator_each_sublatticeの対称性管理のためのモデル構造に使われるので多元系でちゃんと確認。
			for index in list(content):
				self.tag2grouped_basis[prim_cell[index].tag] = group

	@property
	def index_by_basis(self):
		return self.atoms_mng.basisindex_by_tag(self.tag2grouped_basis)

	@property
	def all_cf_names(self):
		num_bf = len(self.basis_functions)
		return self.cluster_list.get_all_cf_names(num_bf)

	@property
	def cluster_mng(self):
		if self._cluster_mng is None:
			kwargs = {}
			if not self.include_background_atoms:
				kwargs["background_syms"] = self.get_bg_syms()
			self._cluster_mng = ClusterManager_each_sublattice(self.prim_cell, **kwargs)
		return self._cluster_mng

	def create_cluster_list_and_trans_matrix(self):
		"""Prepares the internal cache objects by calculating cluster related properties"""
		self._cluster_list, self._trans_matrix = self.cluster_mng.build_all(
			self.atoms,
			self.max_cluster_dia,
			self.index_by_sublattice,
		)

#import output_symmetry
import numpy as np
from cython_process import C_position_to_index 

def cluster2index(config,prim,clusters,target_structure): ##おそらくtarget_structureをおおきくすれば大きなQUBOも作れる(周期境界を超えるかどうかをtarget_structure基準で見ていそう)
	rotations,translations = output_symmetry(config,None,target_structure)
	symmetry_equivalent_index = symmetry_equivalent_arrange(target_structure,{"rotations":rotations,"translations":translations}).T
	size = np.sqrt(np.sum(target_structure.cell**2,axis=1))/np.sqrt(np.sum(prim.cell**2,axis=1))
	cluster_index_dict = {}
	for cluster in clusters:
		if cluster.name == "c0":
			cluster_index_dict["c0"] = ""
		elif "c1" in cluster.name: ##もっと大きなクラスターも取れるが、indicesに入っていたのでどうでもよくなった
			index_list = []
			for figure in cluster.figures:
				figure_positions = []
				for fourvector in figure.components:
					figure_positions.append((fourvector.to_scaled(prim)/size).tolist())
				figure_positions.append([0,0,0])
				figure_positions = np.array(figure_positions)
				figure_positions = (figure_positions-np.min(figure_positions,axis=0)//1)[:-1]
				figure_positions = np.where(figure_positions==1.0,0.0,figure_positions)
				figure_index = C_position_to_index(figure_positions,target_structure.get_scaled_positions()[np.newaxis,:,:])
				index_list.append(figure_index[:,0].tolist())
			cluster_index_dict[cluster.name] = np.unique(symmetry_equivalent_index[figure_index[0,0]]).reshape(-1,1).tolist()
		else:
			cluster_index_ = np.array(cluster.indices)
			cluster_index  = np.sort(cluster_index_, axis=1)
			all_index_     = symmetry_equivalent_index[np.unique(cluster_index,axis=0).tolist()].transpose((1,0,2)).reshape(cluster.size,-1)
			all_index      = np.sort(all_index_, axis=0)
			cluster_index_dict[cluster.name] = np.unique(all_index,axis=1).T.tolist()
	return cluster_index_dict

def cf_calculator_each_sublattice(config,atoms,settings,cf_dict,cluster_list):
	model = atoms.copy()
	model.set_atomic_numbers(model.get_tags())
	cluster_index_dict = cluster2index(config,settings.cluster_mng.prim,cluster_list.clusters,model) # こっちだとxy方向上に別のサイトがいる場合にも対応できるが、bulkのようにsizeより小さいモデルが出てきたときにややこしい->そもそもindexが違う
	basis              = settings.basis_functions ##現時点では1種類の吸着サイトについて考えている。複数あると複数の基底のproductをやる必要がある。平均の取り方がわからないので一旦放置
	substitute_ratio   = config('preparation','substitute_ratio',list)
	if len(substitute_ratio)!=1:
		raise NotImplementedError("This version of SurfQit does not support the cluster expansion for multi sites.")
		# for one site, basis is like [{'Au': 1.2649110640673518, 'Cu': -1.2649110640673518, 'Pd': 0.6324555320336759, 'Pt': -0.6324555320336759}, {'Au': 1.0, 'Cu': 1.0, 'Pd': -1.0, 'Pt': -1.0}, {'Au': 0.6324555320336764, 'Cu': -0.6324555320336764, 'Pd': -1.2649110640673515, 'Pt': 1.2649110640673515}]
		# 全体のgrouped basisを見ながらやる必要がある。多分ね
	for cluster_name in cf_dict.keys():
		if cluster_name == "c0": # もとから0が入っている
			continue
		else:
			indices = []
			for cluster_list_name in cluster_list.names:
				if cluster_list_name in cluster_name:
					indices = indices + cluster_index_dict[cluster_list_name]
			indices = np.unique(indices,axis=0)
			cf = 0
			for index in indices:
				each_cf = 1.0
				for basis_index,atom in zip(list(cluster_name.split("_")[-1]),atoms[index].symbols):
					each_cf *= basis[int(basis_index)][atom]
				cf += each_cf
			cf_dict[cluster_name] = float(cf)/len(indices)
	return cf_dict

import glob
import numpy as np
import ase.io
# import cf_calculator_each_sublattice

def validate_calculatedCF(cluster_expansion_excutor): ##没関数, each layerではcf_matrixとこの中身が同じ機序で動いているはずだが、意外と検証で違う値が出てきたりする。evaluatorの中は自作ではないため
	ce                 = cluster_expansion_excutor
	reference_matrix   = ce.evaluator.cf_matrix
	CF_matirx          = np.zeros_like(reference_matrix)
	all_path           = glob.glob(ce.CE_traindata_path+"id*.cif")
	id_list            = []
	for path in all_path:
		if "afteropt" not in path:
			id_list.append(int(path[path.find("id")+2:-4]))
	id_list.sort()
	for index,model_id in enumerate(id_list):
		model = ase.io.read(ce.CE_traindata_path+f"id{model_id}.cif")
		print("In this function (validate_calculatedCF), tags of loaded model will assign along z height.")
		print("If you require other behavior, please report us on GitHub.")
		height_list = np.unique(model.positions[:,2].round(decimals=5))
		tag         = model.positions[:,2].round(decimals=5)
		for i,height in enumerate(height_list):
			tag     = np.where(tag==height,i,tag)
		model.set_tags(list(map(int,tag)))
		cf_dict = {name: 1.0 for name in ce.evaluator.cf_names}
		cf_dict = cf_calculator_each_sublattice(ce.config,model,ce.settings,cf_dict,ce.settings.cluster_list)
		CF_matirx[index] = np.array([value for value in cf_dict.values()])
	return reference_matrix,CF_matirx

from clease.cluster import ClusterManager
from itertools import product
from clease.cluster import Cluster
from clease.datastructures import Figure, FourVector
from clease.cluster.cluster_fingerprint import ClusterFingerprint

class ClusterManager_each_sublattice(ClusterManager):
	def __init__(self, prim_cell, background_syms = None):
		super().__init__(prim_cell, background_syms)

	def build(self, max_cluster_dia):
		if not self.requires_build(max_cluster_dia):
			return
		self._prepare_new_build(max_cluster_dia)

		num_lattices = range(len(self.prim))
		all_fps = []
		all_figures = []
		lattices = []

		for ref_lattice, (indx, diameter) in product(num_lattices, enumerate(max_cluster_dia)):
			cluster_size = indx + 2
			figures, fps = self.generator.generate(cluster_size, diameter, ref_lattice)

			all_fps += fps
			all_figures += figures
			lattices += [ref_lattice] * len(figures)

		names = self._get_names(all_fps)
		for figures, fp, name, ref_lattice in zip(all_figures, all_fps, names, lattices):
			cluster_size = figures[0].size
			diameter = figures[0].get_diameter(self.prim)
			eq_sites = self.generator.equivalent_sites(figures[0])
			
			name_index_relation = {}
			for figure in figures: ##sublattice毎にfigureの名前を分ける。forが多いので重たくなりそう。あまりに遅いなら考える
				##sublatticeとは結局sizeで広げる前のモデルのindexなので、モデル中に等価なサイトが入っていないと仮定している。
				##完全等価なら従来のClusterManagerを使えばよい。ところどころ等価ならあきらめて分けるか、tagで分けるか。(self.generator.generateをいじればいいが未実装)
				sublattice_list = []
				for fourvector in figure.components:
					sublattice_list.append(str(fourvector.sublattice))
				sublattice_list.sort()
				if "_".join(sublattice_list) not in name_index_relation.keys():
					name_index_relation["_".join(sublattice_list)] = [figure]
				else:
					name_index_relation["_".join(sublattice_list)].append(figure)
			
			for key,value in name_index_relation.items():
				cluster = Cluster(
					name=name+"__"+key,
					size=cluster_size,
					diameter=diameter,
					fingerprint=fp,
					figures=value,
					equiv_sites=eq_sites,
					group=ref_lattice, # ここのgroupはprim内のatom毎に割り振られている(等価かどうかは関係がない)
				)
				self.clusters.append(cluster)

		for index,atom in enumerate(self.prim):
			self.clusters.append(
				Cluster(
					name="c1"+"__"+str(atom.tag),
					size=1,
					diameter=0.0,
					fingerprint=ClusterFingerprint([1.0]),
					figures=[Figure([FourVector(0, 0, 0, index)])],
					equiv_sites=[],
					group=index, # groupはprimのindexについて振られている。等価なサイトを考えるためにnameはtagで判別した。Figureは全てのサイトを基底関数に含めるためにむしろtagではなくindexでみる必要がある
				)
			)

		self.clusters.append(
			Cluster(
				name="c0",
				size=0,
				diameter=0.0,
				fingerprint=ClusterFingerprint([0.0]),
				figures=[],
				equiv_sites=[],
				group=0,
			)
		)

		self.clusters.sort()

	def cluster_list_for_template(self, template):
		unique = self.unique_four_vectors()
		lut = self.fourvec_to_indx(template, unique)
		ref_indices = [lut[FourVector(0, 0, 0, i)] for i in range(self.generator.num_sub_lattices)]

		template_cluster_list = deepcopy(self.clusters)
		for cluster in template_cluster_list:
			cluster.ref_indx = int(ref_indices[cluster.group])
			if cluster.size in {0, 1}:
				cluster.indices = []
			else:
				cluster.indices = self.generator.to_atom_index(cluster, lut)
		return template_cluster_list

import numpy as  np
from clease.settings.concentration import Concentration
#import bulk2slab
#import adjust_concentration_each_layer
#import ClusterExpansionSettings_nowarning

def CESlab_nowarning(conventional_cell,miller,basis_elements,ranges,grouped_basis,vacuum_thickness,model_type,size,**kwargs):
	"""##変える, いろいろ変えたので、最初の方は動作確認が必要。意図したfigureが取れているかとか。そういう関数を作る
	:param conventional_cell:
		Bulk lattice structure. Note that the unit-cell
		must be the conventional cell - not the primitive cell. One can also
		give the chemical symbol as a string, in which case the correct bulk
		lattice will be generated automatically.

	:param miller:
		Surface normal in Miller indices (h,k,l).

	:param concentration:
		Class for restricting the concentrations

	:param size:
		Size of the simulations cell. The third number represents the number of
		layers. The two first are repetitions of the in-plane unit vectors

	For more kwargs, see docstring of :class:`clease.settings.ClusterExpansionSettings`.
	"""
	if model_type == "slab_from_bulk":
		prim,size = bulk2slab(conventional_cell,miller,size,vacuum_thickness)
		if grouped_basis[0] == "all":
			grouped_basis = [tuple(range(len(prim)))] # 前はこうなっていた[tuple(range(len(conventional_cell)))]
		basis_elements_each_layer,grouped_basis_each_layer,ranges_each_layer = adjust_concentration_each_layer(prim,basis_elements,grouped_basis,ranges)
		concentration = Concentration(basis_elements=basis_elements_each_layer,grouped_basis=grouped_basis_each_layer)
		concentration.set_conc_ranges(ranges=ranges_each_layer)
	else:
		prim = conventional_cell
		raise NotImplementedError("動作未検証, 酸化物など、ASEで切るのがめんどくさいやつ用だが、grouped_basisを設定したほうが良い。")
		if grouped_basis[0] == "all":
			grouped_basis = [tuple(range(len(conventional_cell)))]
		concentration = Concentration(basis_elements=basis_elements,grouped_basis=grouped_basis)
		concentration.set_conc_ranges(ranges=ranges)

	# Slab should always have one cell vector along the z-axis
	settings = ClusterExpansionSettings_nowarning(prim, concentration, size=size, **kwargs)

	dict_rep = conventional_cell.todict()
	for k, v in dict_rep.items():
		if isinstance(v, np.ndarray):
			dict_rep[k] = v.tolist()

	settings.kwargs.update(
		{"factory": "CESlab_nowarning","miller": miller,"conventional_cell": dict_rep,"size": size}
	)
	return settings

from ase.db import connect
from ase.calculators.emt import EMT
import pandas as pd
import numpy as np
import ase.io
from ase.optimize import LBFGSLineSearch
from clease.tools import update_db
#import CSV_Reader

class Calculator():
	def __init__(self,config):
		self.storage_path = config('preparation','storage_path',str)
		self.db_name      = self.storage_path+"train_data/training_data.db"
		self.csv_path     = self.storage_path+"train_data/energy_data.csv"
		self.method       = config('calculation','calculate_method',list)
		self.fmax         = config('calculation','fmax',float)
		if self.method[0] == "EMT":
			self.calc = EMT()
			self.reader = CSV_Reader(config)
		elif self.method[0] == "csv_reader":
			self.reader = CSV_Reader(config)
		else:
			raise ValueError(f"{self.method[0]}is not implemented")

	def __call__(self,data,calculate_function=None, constraint_dict=None):
		model_id = data.id
		model = data.toatoms()
		calculate_function(model,model_id,constraint_dict)
		
	def optimization(self,model,model_id,constraint_dict=None): ##複雑な計算をやりたいときは外部で計算してcsvを読んだ方が良い
		model.calc = self.calc
		if constraint_dict is not None:
			model.set_constraint(constraint_dict["constraint"])
			if "filter" in constraint_dict:
				model = constraint_dict["filter"]["function"](model,*constraint_dict["filter"]["args"])
		opt = LBFGSLineSearch(model, logfile=None)
		opt.run(fmax=self.fmax)
		ase.io.write(self.storage_path+f"train_data/id{str(model_id)}_afteropt.cif",model,format="cif")
		model = ase.io.read(self.storage_path+f"train_data/id{str(model_id)}_afteropt.cif",format="cif") ##calculatorによってはopt後と再読み込み後でのget_potential_enegyが違うことがあるので再現性のある値を取りたい。
		self.single_point(model,model_id)

	def single_point(self,model,model_id,constraint_dict=None):
		model.calc = self.calc
		energy = model.get_potential_energy()
		self.csv_write(model,model_id,energy)

	def csv_write(self,model,model_id,energy,constraint_dict=None):
		energy_per_atom = energy/len(model)
		if len(np.unique(model.symbols))>1:
			elements = np.unique(model.symbols)
			atom_numbers = np.zeros_like(elements,dtype=np.int32)
			for atom in model:
				atom_numbers[np.where(elements==atom.symbol)[0]] += 1
		else:
			elements = np.unique(model.symbols)
			atom_numbers = [len(model)]
		df = pd.read_csv(self.csv_path)
		column = df.columns
		additional_data = np.zeros(len(column))
		additional_data[0] = model_id
		additional_data[1] = energy_per_atom
		additional_data[2] = np.nan
		for element_index,element in enumerate(elements):
			additional_data[column.get_loc(element)] = atom_numbers[element_index] 
		df = pd.DataFrame(additional_data.reshape(1,-1),index=[model_id])
		df.to_csv(self.csv_path,mode="a",header=False,index=False)

	def csv_read(self,data):
		model_id = data.id
		model = data.toatoms()
		self.reader.set_model_id(model_id)
		model.calc = self.reader
		model.get_potential_energy()
		update_db(uid_initial=model_id, final_struct=model, db_name=self.db_name)

	def csv_write_formationE(self):
		df = pd.read_csv(self.csv_path)
		reference_energies = {}
		for energy,composition in zip(df[df.columns[1]],df[df.columns[3:]].to_numpy()):
			nonzero_index = np.where(composition!=0)[0]
			if len(nonzero_index) == 1:
				reference_energies[df.columns[3:][nonzero_index].to_numpy()[0]]=energy
		for row in df.itertuples():
			formation_energy = row.energy_per_atom * sum(row[4:])
			for key,value in reference_energies.items():
				formation_energy -= row[df.columns.to_list().index(key)+1] * value #/ sum(row[3+1:]) ## Cleaseの内部で原子数で割り付けていそう
			df.at[row.Index,"formation_energy"] = formation_energy
		df.to_csv(self.csv_path,mode="w",header=True,index=False)

#optのconstraintはどのようにするか考える*args?
##subprocess.run(list(.split(" ")))←動くか富岳でチェック, shell=Trueもできなくはないけどセキュリティリスク

from ase.calculators.calculator import Calculator as ASE_Calculator
import pandas as pd

all_changes = ['positions', 'numbers', 'cell', 'pbc',
               'initial_charges', 'initial_magmoms']

class CSV_Reader(ASE_Calculator): # based on EMT
    implemented_properties = ['energy']
    nolabel = True
    default_parameters = {'asap_cutoff': False}

    def __init__(self, config, **kwargs):
        storage_path = config('preparation','storage_path',str)
        self.csv_path     = storage_path+"train_data/energy_data.csv"
        ASE_Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        self.df = pd.read_csv(self.csv_path)
        ASE_Calculator.calculate(self, atoms, properties, system_changes)
        self.results['energy'] = float(self.df[self.df["id"]==self.id]["formation_energy"].values)

    def get_potential_energy(self,atoms):
        return self.get_property("energy", atoms)

    def set_model_id(self,model_id):
        self.id = model_id

from clease.structgen import NewStructures
# import CorrFunction_each_sublattice

class NewStructures_each_sublattice(NewStructures):
	def __init__(self,config,settings,generation_number = None,struct_per_gen = 5,check_db = True):
		self.check_db = check_db
		self.settings = settings
		self.config   = config
		self.corrfunc = CorrFunction_each_sublattice(self.config,self.settings)
		self.struct_per_gen = struct_per_gen

		if generation_number is None:
			self.gen = self._determine_gen_number()
		else:
			self.gen = generation_number

	def generate_probe_structure(self, atoms=None, init_temp=None, final_temp=None, num_temp=5, num_steps_per_temp=1000, approx_mean_var=True, num_samples_var=10000):
		if not approx_mean_var:
			if not os.path.isfile("probe_structure-sigma_mu.npz"):
				self._generate_sigma_mu(num_samples_var)

		current_count = 0
		num_attempt = 0
		num_to_generate = self.num_to_gen()

		while current_count < num_to_generate:
			self.settings.set_active_template(atoms=atoms)
			num_struct = self.num_in_gen()
			if num_struct >= self.struct_per_gen:
				break

			struct = self._get_struct_at_conc(conc_type="random")

			ps = ProbeStructure_each_sublattice(self.settings,self.config,struct,init_temp,final_temp,num_temp,num_steps_per_temp,approx_mean_var)
			probe_struct, cf = ps.generate()

			probe_struct.calc.results.pop("energy")
			formula_unit = self._get_formula_unit(probe_struct)
			if self._exists_in_db(probe_struct, formula_unit):
				num_attempt += 1
				if num_attempt >= max_attempt:
					raise MaxAttemptReachedError(msg)
			else:
				num_attempt = 0

			self.insert_structure(init_struct=probe_struct, cf=cf)
			current_count += 1

from clease.structgen import ProbeStructure
# import CorrFunction_each_sublattice

class ProbeStructure_each_sublattice(ProbeStructure):
	def __init__(self,config,settings,atoms,init_temp=None,final_temp=None,num_temp=5,num_steps_per_temp=1000,approx_mean_var=False):
		ProbeStructure.__init__(self, settings, atoms, init_temp, final_temp, num_temp, num_steps_per_temp, approx_mean_var)
		self.corrFunc = CorrFunction_each_sublattice(config,self.settings)

from clease.corr_func import CorrFunction
from ase import Atoms
# import cf_calculator_each_sublattice

class CorrFunction_each_sublattice(CorrFunction):
	def __init__(self,config,settings):
		self.config   = config
		self.settings = settings

	def get_cf_by_names(self, atoms, cf_names):
		if isinstance(atoms, Atoms):
			self.set_template(atoms)
		else:
			raise TypeError("atoms must be Atoms object")

		self._confirm_cf_names_exists(cf_names)

		cf = {name: 1.0 for name in cf_names}
		cf = cf_calculator_each_sublattice(self.config,atoms,self.settings,cf,self.settings.cluster_list)
		return cf

def check_position_overlap(new,reference,tol = 1e-1):
	if np.max(np.abs(new.positions-reference.positions)) > tol:
		raise ValueError("The input model is not the same structure of the target model.")

import numpy as  np

def check_max_cluster_size(prim_model,size,max_cluster_size):
	adjacent_mask = np.zeros((3,3,3,3))
	adjacent_mask[0,:,:,0] = -1
	adjacent_mask[2,:,:,0] = 1
	adjacent_mask[:,0,:,1] = -1
	adjacent_mask[:,2,:,1] = 1
	adjacent_mask[:,:,0,2] = -1
	adjacent_mask[:,:,2,2] = 1
	model = prim_model.repeat(size)
	adjacent_length = np.sum((adjacent_mask*np.sqrt(np.sum(model.cell**2,axis=1))[np.newaxis,np.newaxis,np.newaxis,:])**2,axis=3)
	adjacent_length[1,1,1] = np.max(adjacent_length)
	adjacent = np.argwhere(adjacent_length==np.min(adjacent_length))[0]-1
	adjacent_translate = adjacent*np.sqrt(np.sum(model.cell**2,axis=1))
	adjacent_translate_length = np.sqrt(np.sum(adjacent_translate**2))
	adjacent_prim = prim_model.positions+adjacent_translate
	extended_prim = prim_model.repeat(size).repeat((3,3,3)).positions
	adjacent_distance = np.sqrt(np.sum((adjacent_prim[:,np.newaxis,:]-extended_prim[np.newaxis,:,:])**2,axis=2)) # ここのまえにあまりに長すぎる距離(真逆とか)は丸めてもいいかも。2乗がえぐそう
	prim_distance     = np.sqrt(np.sum((prim_model.positions[:,np.newaxis,:]-extended_prim[np.newaxis,:,:])**2,axis=2))
	distance          = np.concatenate([adjacent_distance[np.newaxis,:,:],prim_distance[np.newaxis,:,:]])
	distance          = np.max(distance,axis=0)
	criterion         = np.min(distance)
	if max(max_cluster_size)>=criterion: ## デバッグ的にはないほうがいいかも
		raise ValueError(f"The max_cluster_size should be up to {criterion} Å (now: {max(max_cluster_size)} Å). This setting leads to calculate ECIs for the same atom configuration between 2 or more clusters due to the periodic boundary condtion. Please reduce max_cluster_size or increase model size.")

import os
import numpy as np
import ase.io 
from ase import Atom
from ase.build.tools import sort

class Slab_Processor:
	def __init__(self,cif_list=None,model_type=None,vacuum=None,config=None,leyer_number=None):
		if cif_list is not None:
			self.cif_list         = cif_list
		self.radiusdict={'Au':1.44,'Co':1.25,'Pd':1.37,'Zn':1.37,'Fe':1.26,'P':1.1,'Ce':1.82,'Pt':1.39,'O':0.74,'Ba':2.24,'Cu':1.28,'Ni':1.25,'Ru':1.34}
		self.atomcandidate = ["H","C","N","O"]
		if model_type is not None:
			self.set_model_type(model_type)
		if vacuum is not None:
			self.set_vaccum_thickness(vacuum)
		if config is not None:
			self.vacuum_thickness = config('preparation','vacuum',float)
		if leyer_number is not None:
			self.set_leyer_number(leyer_number)

	def __call__(self,function_list,profix=None, fileformat="cif"):
		self.load()
		for function in function_list:
			self.process_iter(function)
		if fileformat == 'cif':
			self.save_cif(profix)
		elif fileformat == 'vasp':
			self.save_vasp(profix)

	def update_radiusdict(self,symbol,radius):
		self.radiusdict[symbol] = radius

	def selective_dynamics(self,POSCAR_list,adsorbate_number,fix_relax_list):
		for POSCAR in POSCAR_list:
			model = ase.io.read(POSCAR)
			split_index = int((len(model)-adsorbate_number) * fix_relax_list[0]/sum(fix_relax_list))
			fixindex = np.argsort(model.positions[:-adnum,2])[:split_index] # 上下吸着物がにいても取ってこれる
			relaxindex = np.argsort(model.positions[:-adnum,2])[split_index:]
			adindex= np.argsort(model.positions[-adnum:,2])
			with open(resultdir+'/POSCAR', "r") as f:
				lines = f.readlines()
			if fix_relax_layer[0]>0:
				lines.insert(7, 'S'+'\n')
				for i in range(9,len(lines)):
					if i in fixindex+9:
						lines[i] = lines[i].rstrip() + ' F F F \n'
					elif i in relaxindex+9:
						lines[i] = lines[i].rstrip() + ' T T T \n'
					elif i in adindex+9:
						lines[i] = lines[i].rstrip() + ' T T T \n'
			with open(resultdir+'/POSCAR', "w") as f:
				f.writelines(lines)

	def set_cif_list(self,cif_list):
		self.cif_list = cif_list

	def set_vaccum_thickness(self,vacuum_thickness):
		self.vacuum_thickness = vacuum_thickness

	def set_leyer_number(self,leyer_number):
		self.leyer_number = leyer_number

	def set_model_type(self,model_type):
		if model_type == "slab":
			self.bottom_bool = False
		elif model_type == "sheet":
			self.bottom_bool = True

	def load(self):
		self.models = []
		for cif in self.cif_list:
			model = ase.io.read(cif)
			self.models.append(model)

	def save_cif(self,profix=None):
		for cif_path,model in zip(self.cif_list,self.models):
			if profix is not None:
				ase.io.write(f"{cif_path[:-4]}_{profix}.cif",model)
			else:
				ase.io.write(f"{cif_path[:-4]}.cif",model)

	def save4visualize(self,profix=None):
		for cif_path,model in zip(self.cif_list,self.models):
			element_dict={}
			for atom in model:
				if atom.symbol not in element_dict.keys():
					element_dict[atom.symbol] = 1
				else:
					element_dict[atom.symbol] += 1
			name = ""
			element_list = list(element_dict.keys())
			element_list.sort()
			for element in element_list:
				name = name+element+str(element_list[element])
			if profix is not None:
				ase.io.write(f"{os.path.dirname(cif_path[:-4])+'/'+name}_{profix}.cif",model)
			else:
				ase.io.write(f"{os.path.dirname(cif_path[:-4])+'/'+name}.cif",model)

	def save_vasp(self,profix=None):
		for cif_path,model in zip(self.cif_list,self.models):
			if profix is not None:
				ase.io.write(f"{cif_path[:-4]}_{profix}.vasp",model)
			else:
				ase.io.write(f"{cif_path[:-4]}.vasp",model)

	def process_iter(self,function):
		for index,model in enumerate(self.models):
			self.models[index] = function(model)

	def sort(self,model):
		model = sort(model)
		return model

	def remove_vacuum(self,model):
		cell = model.get_cell()
		cell[2,:] = cell[2,:]*(cell[2,2]-self.vacuum_thickness)/cell[2,2]
		model.set_cell(cell)
		return model

	def add_vacuum(self,model):
		cell = model.get_cell()
		cell[2,:] = cell[2,:]*(cell[2,2]+self.vacuum_thickness)/cell[2,2]
		model.set_cell(cell)
		return model

	def undo_z_wrap(self,model):
		z = model.get_cell()[2,2]
		for atom in model:
			if atom.position[2] > z*(self.leyer_number-0.9)/self.leyer_number:
				atom.position[2] -= z
		offset = np.min(model.positions[:,2])
		for atom in model:
			atom.position[2] -= offset
		return model

	def adsorption_Hydrogen(self,model):
		return self.adsorption_moledule(model,[[0,0,0,2.1]])

	def adsorption_Oxygen(self,model):
		return self.adsorption_moledule(model,[[3,0,0,2.1]])

	def adsorption_Carbon(self,model):
		return self.adsorption_moledule(model,[[1,0,0,2.1]])

	def adsorption_CO(self,model):
		return self.adsorption_moledule(model,[[1,0,0,2.1],[3,0,0,2.1+1.18]])

	def adsorption_moledule(self,model,adsorbate_data):
		surface_judge = self.search_surfaceatoms(model)
		top, bottom = self.get_surfaceatoms(model, surface_judge)
		admodel= model.copy()
		vacuum = model.copy()
		del vacuum[:]
		count=0
		for i in top:
			for atom in adsorbate_data:
				vacuum+=Atom(self.atomcandidate[atom[0]],model[i].position+np.array(atom[1:]))
				count+=1
		if self.bottom_bool:
			for i in bottom:
				for atom in adsorbate_data:
					vacuum+=Atom(self.atomcandidate[atom[0]],model[i].position-np.array(atom[1:]))
					count+=1
		elementlist=np.unique(vacuum.get_chemical_symbols())
		for element in elementlist:
			for i in vacuum:
				if i.symbol==element:
					admodel+=Atom(i.symbol,i.position)
		return admodel

	def search_surfaceatoms(self,model):
		model, index, original_indices = self.model_expansion(model)
		atomradiuslist=[self.radiusdict[i.symbol] for i in model]
		atompositionlist=[i.position for i in model]
		original_zlist=np.arange(-1,1.2,0.2)
		surfacelist=[]
		for n,i in enumerate(model):
			target_radius=self.radiusdict[i.symbol]+1.11
			zlists=np.repeat(target_radius*original_zlist,16).reshape(len(original_zlist)*16,1)
			xlists=np.cos(np.linspace(0, 2.0*np.pi, 16))
			ylists=np.sin(np.linspace(0, 2.0*np.pi, 16))
			xylists=np.array(list(np.concatenate((xlists.reshape(len(xlists),1),ylists.reshape(len(ylists),1)),axis=1))*len(original_zlist))
			ratiolist=(target_radius)**2-(zlists)**2
			balllist=np.hstack((np.sqrt(ratiolist)*xylists,zlists))

			positions=i.position+balllist
			repeat_positions=np.repeat(positions.reshape(len(positions),1,3),len(atompositionlist),axis=1)
			distance=np.sqrt(np.sum((repeat_positions-atompositionlist)**2,axis=2))
			totalradius=np.array(atomradiuslist)+1.1
			if True in np.all(distance > totalradius,axis=1):
				surfacelist.append(n)

		surface_judge=np.zeros(len(model))
		for i in original_indices:
			if i in surfacelist:
				surface_judge[index[i]]=1
		return surface_judge

	def model_expansion(self,model):
		expand_ratio=(3,3,1)
		index = np.tile(np.arange(len(model)),np.prod(expand_ratio))
		originalnum = int(np.prod(expand_ratio) / 2)
		original_indices = np.arange(originalnum * len(model), (originalnum + 1) * len(model))
		model = model * expand_ratio
		return model, index, original_indices

	def get_surfaceatoms(self,model,connectivity):
		surf_atoms = connectivity != 0
		indices = np.argwhere(surf_atoms).flatten()
		relative_zcoords= model.positions[indices,2] - np.mean(model.positions[:,2])
		top = indices[relative_zcoords >= 0]
		bottom = indices[relative_zcoords < 0]
		return top, bottom

import numpy as np
import pandas as pd
import itertools
import random
import glob
from amplify_bbopt import blackbox, BinaryVariableList, DatasetGenerator, KernelQAOptimizer
import amplify_bbopt # equal_toのため。フォルダ分けでアニーリングのequal_toと被らなくなれば上に統合できる
import re
import ase.io
from ase import Atoms
from ase.db import connect
from ase.build import add_vacuum
from scipy.special import comb
import json
import os
from clease.settings.concentration import Concentration
from clease.settings import CECrystal,CESlab
from clease.basis_function import BasisFunction, BinaryLinear, Polynomial, Trigonometric
#from clease.settings.settings_slab import remove_vacuum_layers
from clease.structgen import NewStructures
from clease import Evaluate
#import CESlab_nowarning
#import Calculator
#import NewStructures_each_sublattice
# import Fixstars_Amplify

##去年の時点で遅いとされていた実装などでもベーシックなら入れているので、今後元を参考にしつつCythonizeしていけ
class ClusterExpansionExcutor():
	def __init__(self,config,seed=None,confirm_max_cluster_size=True):
		"""
		input:
			config: Prepare_Configクラスのもの。
			seed: 構造生成をそろえて比較したいとき用
		output:
			None
		description:
			CEの準備をします。
		"""
		##こんなにいるか?フォルダ構成は考える, bulkはチェックが必要
		self.config               = config
		self.storage_path         = config('preparation','storage_path',str)
		self.db_name              = config('preparation','storage_path',str)+"train_data/training_data.db"
		self.host                 = ase.io.read(self.storage_path+'host.cif',format='cif')
		self.CE_traindata_path    = self.storage_path+"train_data/"
		self.model_type           = config('preparation','model_type',str)
		self.miller_index         = config('preparation','miller',tuple)
		self.vacuum_thickness     = config('preparation','vacuum',float)
		self.concentration        = config('preparation','substitute_ratio',list)
		self.grouped_basis        = config('preparation','group_index',list)
		self.supercell_multiple   = config('preparation','supercell_multiple',list)

		self.basis_function       = config('cluster','basis_function',str)
		self.max_cluster_size     = config('cluster','max_cluster_size',list)
		self.basis_elements       = []

		self.parallel_calculation = config('calculation','train_data_num',int)
		self.calculation_method   = config('calculation','calculate_method',list)
		self.scoring_scheme       = config('calculation','scoring_scheme',str)
		self.nsplits              = config('calculation','split_number',int)
		self.fitting_scheme       = config('calculation','fitting_scheme',str)
		self.alpha                = config('calculation','alpha',list)
		self.annealing_machine    = config('calculation','annealing_machine',str)

		self.confirm_cluster_size = confirm_max_cluster_size
		self.calculator           = Calculator(self.config)
		if seed is not None:
			random.seed(int(seed))
			np.random.seed(int(seed))
		os.makedirs(self.CE_traindata_path,exist_ok=True)
		ranges         = []
		for substitution in self.concentration:
			self.basis_elements.append(list(substitution.keys()))
			ranges.append(list(map(tuple,substitution.values())))
		if len(self.supercell_multiple)==3:
			if self.model_type == "bulk":
				if self.grouped_basis[0] == "all":
					grouped_basis = [tuple(range(len(self.concentration)))]
				concentration = Concentration(basis_elements=self.basis_elements,grouped_basis=grouped_basis)
				concentration.set_conc_ranges(ranges=ranges)
				output_symmetry(self.config,host=self.host)
				with open(self.storage_path+'symmetry.json','r') as f:
					symmetry_data = json.load(f)
				self.settings = CECrystal(concentration,spacegroup=str(symmetry_data['space group symbol']),\
											basis=self.host.get_scaled_positions(),cell=self.host.get_cell()[:],\
											basis_func_type=self.basis_function,db_name=self.db_name,\
											size=self.supercell_multiple,max_cluster_dia=self.max_cluster_size,\
											crystal_kwargs={"onduplicates": "keep"}) # no crystal_kwargs leads to print user warning {symmetry sites in self.host} times. 
				if confirm_max_cluster_size:
					check_max_cluster_size(self.settings.cluster_mng.prim,self.supercell_multiple,self.max_cluster_size)
			else:
				if self.vacuum_thickness<max(self.max_cluster_size): ## slab from bulkのみしかまだ考えていない
					print("Warning: vacuum @ config should be larger than maximum max_cluster_size @ config to avoid generating interaction between topmost layer and bottom layer through the vacuum layer.")
					print("vacuum @ config set to maximum max_cluster_size+1 @ config.")
					self.vacuum_thickness = max(self.max_cluster_size)+1
					##config自体も書き換えたほうがユーザーライク?
				self.settings = CESlab_nowarning(self.host,self.miller_index,self.basis_elements,ranges,self.grouped_basis,vacuum_thickness=self.vacuum_thickness,model_type=self.model_type,size=(self.supercell_multiple[0],self.supercell_multiple[1],self.supercell_multiple[2]),\
											basis_func_type=self.basis_function,db_name=self.db_name,max_cluster_dia=self.max_cluster_size)
				if confirm_max_cluster_size:
					check_max_cluster_size(self.settings.cluster_mng.prim,(self.supercell_multiple[0],self.supercell_multiple[1],1),self.max_cluster_size)
		elif len(self.supercell_multiple)==1:
			if self.model_type == "bulk":
				if self.grouped_basis[0] == "all":
					grouped_basis = [tuple(range(len(self.concentration)))]
				concentration = Concentration(basis_elements=self.basis_elements,grouped_basis=grouped_basis)
				concentration.set_conc_ranges(ranges=ranges)
				output_symmetry(self.config,host=self.host)
				with open(self.storage_path+'symmetry.json','r') as f:
					symmetry_data = json.load(f)
				self.settings = CECrystal(concentration,spacegroup=str(symmetry_data['space group symbol']),\
											basis=self.host.get_scaled_positions(),cell=self.host.get_cell()[:],\
											basis_func_type=self.basis_function,db_name=self.db_name,\
											supercell_factor=self.supercell_multiple[0],max_cluster_dia=self.max_cluster_size,\
											crystal_kwargs={"onduplicates": "keep"}) # no crystal_kwargs leads to print user warning {symmetry sites in self.host} times. 
				if confirm_max_cluster_size:
					check_max_cluster_size(self.settings.cluster_mng.prim,(self.supercell_multiple[0],self.supercell_multiple[0],self.supercell_multiple[0]),self.max_cluster_size)
			else:
				raise ValueError(f"do not support {len(self.supercell_multiple)}-sized supercell_multiple for {self.model_type} mode")
		else:
			raise ValueError(f"do not support {len(self.supercell_multiple)}-sized supercell_multiple")
		element_candidates = list(np.unique(np.array(sum(self.basis_elements,[]))))
		if not os.path.isfile(self.CE_traindata_path+"energy_data.csv"):
			df = pd.DataFrame(columns = ["id","energy_per_atom","formation_energy"]+element_candidates)
			df.to_csv(self.CE_traindata_path+"energy_data.csv",index=False)

	def reload_config(self,config,cf=None):
		if self.model_type == "bulk":
			raise NotImplementedError("This function with bulk has not implemented yet.")
		symmetry_tolerance = config('preparation','symmetry_tolerance',float)
		if 1e-4 > symmetry_tolerance:
			print("The decimals of cif data is usually 5. The symmetry_tolerance should be a order less than the decimals.")
		os.remove(self.db_name)
		self.__init__(config,confirm_max_cluster_size=self.confirm_cluster_size)
		model_id_list = []
		for candidate in glob.glob(self.CE_traindata_path+"id*.cif"):
			model_id_list.append(re.findall(r'\d+', candidate)[0])
		model_id_list = list(map(int,set(model_id_list)))
		model_id_list.sort()
		for model_id in model_id_list:
			self.add_user_structure(self.CE_traindata_path+f"id{str(model_id)}.cif")
		self.update_formation_energy()

	def model_addition_random(self,generation_number=None):
		"""以前はsublatticeの有無で判別していたが、実際には新しく書いたcf calculatorを使いたいので同じものを使うことにした
		if self.settings.concentration.grouped_basis is None:
			training_data_generator = NewStructures_each_sublattice(self.config,self.settings, generation_number=generation_number, struct_per_gen=self.parallel_calculation)
		elif len(self.settings.concentration.grouped_basis) == 1: ##sublatticeが無ければcfでsublatticeの考慮不要
			training_data_generator = NewStructures(self.settings, generation_number=generation_number, struct_per_gen=self.parallel_calculation)
		else:
			training_data_generator = NewStructures_each_sublattice(self.config,self.settings, generation_number=generation_number, struct_per_gen=self.parallel_calculation)
		"""
		training_data_generator = NewStructures_each_sublattice(self.config,self.settings, generation_number=generation_number, struct_per_gen=self.parallel_calculation)
		training_data_generator.generate_initial_pool()
		db = connect(self.db_name)
		for data in db.select(converged=False):
			atom = data.toatoms()
			ase.io.write(self.CE_traindata_path+f"id{str(data['id'])}.cif",atom,format="cif")

	def add_user_structure(self,model,cf=None):
		if not isinstance(model, Atoms):
			model = ase.io.read(model)
		z = []
		for atom in model:
			if np.round(atom.position[2],decimals=6) not in z:
				z.append(np.round(atom.position[2],decimals=6))
		z.sort()
		for atom in model: # そもそも複数サイトで働かない?⇒basis groupを考える必要。, 酸化物表面も重要。
			atom.tag = z.index(np.round(atom.position[2],decimals=6)) # 高さでtag
		model = wrap_and_sort_by_position(model)
		template = self.settings.template_atoms.weighted_random_template()
		template = wrap_and_sort_by_position(template)
		indices = np.argmin(np.sum((model.positions[np.newaxis,:,:]-template.positions[:,np.newaxis,:])**2,axis=2),axis=1)
		sorted_model = model.copy()[indices]
		check_position_overlap(sorted_model,template)
		sorted_model.positions = template.positions
		if len(self.settings.concentration.grouped_basis) == 1: ##sublatticeが無ければcfでsublatticeの考慮不要
			data_inserter = NewStructures(self.settings, None, struct_per_gen=self.parallel_calculation)
		else:
			data_inserter = NewStructures_each_sublattice(self.config,self.settings, None, struct_per_gen=self.parallel_calculation)
		data_inserter.insert_structure(sorted_model,cf=cf)
		db = connect(self.db_name)
		for data in db.select(converged=False):
			atom = data.toatoms()
			ase.io.write(self.CE_traindata_path+f"id{str(data['id'])}.cif",atom,format="cif")

	def add_further_model(self,struct_per_gen=1,generation_number=None,mode=None,\
						init_temp=2000,final_temp=1, num_temp=10,num_steps_per_temp=5000, random_composition=True, ECI_file_name="ECI.json"):
		training_data_generator = NewStructures_each_sublattice(self.config,self.settings, generation_number=generation_number, struct_per_gen=struct_per_gen)
		if mode == "ground_state":
			template = self.settings.template_atoms.weighted_random_template()
			template = wrap_and_sort_by_position(template)
			eci_dict = {}
			for key,value in zip(*self.get_ECI(ECI_file_name)):
				eci_dict[key] = value
			training_data_generator.generate_gs_structure(atoms=template, init_temp=init_temp, final_temp=final_temp, num_temp=num_temp, num_steps_per_temp=num_steps_per_temp, eci=eci_dict, random_composition=random_composition)
		elif mode == "different":
			raise NotImplementedError("The validation of structure (whether the most far from ex-models or not) has not confirmed yet.")
			training_data_generator.generate_probe_structure()
		else:
			raise ValueError(f"The mode {mode} is not supported in this version of SurfQit.")

	def add_user_structure_origin(self,model): # 今使っていないかも?
		if not isinstance(model, Atoms):
			model = ase.io.read(model)
		model.wrap()
		template = self.settings.template_atoms.weighted_random_template()
		indices = np.argmin(np.sum((model.positions[np.newaxis,:,:]-template.positions[:,np.newaxis,:])**2,axis=2),axis=1)
		sorted_model = model.copy()[indices]
		check_position_overlap(sorted_model,template)
		sorted_model.positions = template.positions
		if len(self.settings.concentration.grouped_basis) == 1: ##sublatticeが無ければcfでsublatticeの考慮不要
			data_inserter = NewStructures(self.settings, None, struct_per_gen=self.parallel_calculation)
		else:
			data_inserter = NewStructures_each_sublattice(self.config,self.settings, None, struct_per_gen=self.parallel_calculation)
		data_inserter.insert_structure(sorted_model)
		db = connect(self.db_name)
		for data in db.select(converged=False):
			atom = data.toatoms()
			ase.io.write(self.CE_traindata_path+f"id{str(data['id'])}.cif",atom,format="cif")
			
	def unconverged_calculate(self,constraint_dict=None):
		db = connect(self.db_name)
		if self.calculation_method[1] != "csv_read":
			for data in db.select(converged=False):
				self.calculate(data,constraint_dict) # 計算を行う
		self.update_formation_energy()

	def update_formation_energy(self):
		self.calculator.csv_write_formationE()
		db = connect(self.db_name)
		for data in db.select(converged=False):
			self.calculator.csv_read(data) # formation_energyでアップデートする

	def calculate(self,data,constraint_dict=None): # return energyをやめた(のでCalculatorもreturn消している) ## 関数説明しろ　dictは{"constraint":[],"filter":{"function":function,"args":[引数たち]}}
		if self.calculation_method[1] == "single_point":
			self.calculator(data,self.calculator.single_point)
		elif self.calculation_method[1] == "optimization":
			self.calculator(data,self.calculator.optimization,constraint_dict)
		else:
			raise ValueError(f"Calculation type {self.calculation_method[1]} is not impremented")

	def ECI_evaluate(self,ECI_file_name="ECI.json",CV_filename=None):
		## polynomial以外はここでエラーが出るはず。reconstructionでは治らなかった
		self.evaluator = Evaluate(settings=self.settings, scoring_scheme=self.scoring_scheme,nsplits=self.nsplits)
		if self.fitting_scheme=="linear":
			self.evaluator.set_fitting_scheme(fitting_scheme=self.fitting_scheme)
		elif CV_filename==None:
			self.evaluator.set_fitting_scheme(fitting_scheme=self.fitting_scheme)
			alpha = self.evaluator.plot_CV(alpha_min=self.alpha[0], alpha_max=self.alpha[1], num_alpha=self.alpha[2])
			self.evaluator.set_fitting_scheme(fitting_scheme=self.fitting_scheme,alpha=alpha)
		else:
			self.evaluator.set_fitting_scheme(fitting_scheme=self.fitting_scheme)
			alpha = self.evaluator.plot_CV(alpha_min=self.alpha[0], alpha_max=self.alpha[1], num_alpha=self.alpha[2],savefig=True,fname=CV_filename)
			self.evaluator.set_fitting_scheme(fitting_scheme=self.fitting_scheme,alpha=alpha)
		self.evaluator.fit()
		cv = self.evaluator.get_cv_score()*1000.0
		if cv >= 5.0:
			print(f"Your CV score is more than 5 meV/atom. ({cv}) You should change CE conditions or/and increase the number of model.")
		self.evaluator.save_eci(self.storage_path+ECI_file_name)
		return cv

	def ECI_evaluate_with_cf_optimization(self, max_cf_num, ECI_file_name="ECI.json",CV_filename=None,init_sample=5,sample=15):
		if self.annealing_machine != "FA": 
			raise NotImplementedError("This version of SurfQit does not support the cf optimization except for Fixstars Amplify.")
		print("CV score by full set of clusters: ", self.ECI_evaluate(ECI_file_name,CV_filename))
		full_cf = np.array(self.get_ECI(ECI_file_name)[0])

		@blackbox
		def objective_lowest_CV(cf_bool_list: list = BinaryVariableList(len(full_cf)))->float:
			self.reload_config(self.config,full_cf[cf_bool_list])
			return self.ECI_evaluate(ECI_file_name,CV_filename)

		objective_lowest_CV.add_constraint(amplify_bbopt.equal_to(objective_lowest_CV.variables.cf_bool_list.to_poly(),float(max_cf_num)))
		generator         = DatasetGenerator(objective=objective_lowest_CV)
		data              = generator.generate(num_samples=init_sample)
		FA                = Fixstars_Amplify(self.config,None,None,None) # QUBOのアニールと設定を共通させるのはいつか不都合が出るかも
		self.cf_optimizer = KernelQAOptimizer(data=data,objective=objective_lowest_CV,client=FA.client)
		self.additional_cf_optimization(ECI_file_name,CV_filename,sample)

	def additional_cf_optimization(self,ECI_file_name,CV_filename,sample):
		full_cf = np.array(self.get_ECI(ECI_file_name)[0])
		self.cf_optimizer.optimize(num_cycles=sample)
		self.reload_config(self.config,full_cf[self.cf_optimizer.best_solution["cf_bool_list"]])
		print("optimized set of clusters: ", full_cf[self.cf_optimizer.best_solution["cf_bool_list"]])
		print("CV score by optimized set of clusters: ", self.ECI_evaluate(ECI_file_name,CV_filename))

	def get_ECI(self,ECI_file_name="ECI.json"):
		self.evaluator.load_eci(self.storage_path+ECI_file_name)
		return self.evaluator.cf_names,self.evaluator.eci

	def get_basis_function(self,unique_element):
		# a basis built for a type of substitution (same element variety)
		if self.basis_function.lower() == "polynomial":
			basis_funcion = Polynomial(unique_element)
		elif self.basis_function.lower() == "trigonometric":
			basis_funcion = Trigonometric(unique_element)
		elif self.basis_function.lower() == "binary_linear":
			basis_function = BinaryLinear(unique_element)
		else:
			msg = f"basis function type {bf_type} is not supported."
			raise ValueError(msg)
		return basis_funcion

	def make_qubit(self,model=None): ## tag設定の観点からここでmodelを受け入れるのは無理?
		if model is None:
			supercell_multiple = self.supercell_multiple
			if self.model_type != "bulk":
				supercell_multiple[2] = 1
			model = self.settings.prim_cell.copy()
			
			model = model.repeat(supercell_multiple)
			model = wrap_and_sort_by_position(model)
			self.QUBO_model = model
		if self.model_type == "bulk":
			print(f"This bulk system has {len(np.unique(model.get_tags()))} types of inequivalent sites, ",*np.unique(model.get_tags()))
			print("If this result is different from your intention, please set tag for each atom of the bulk.")
		elif model.get_tags().sum() == 0.0:
			raise ValueError("The input model does not have specific tag corresponding to symmetrical equivalent sites.")
		if self.grouped_basis == ["all"]:
			grouped_tag = [tuple(set([atom.tag for atom in self.settings.prim_cell]))]
		else: ## まだ直していない
			raise NotImplementedError("This version of SurfQit does not support the multi site.")
			grouped_basis = CE.grouped_basis ## CE.grouped_basisがそもそもtagになっていないが、広げた後に対応させるためにtagに対応させる必要がある。
		basis_dict ={}
		"""ここをtagで管理する""" ## 実際にはCE.basis_elementsは"all"以外でtagに対応してないので外部で処理の必要
		for elements,tags in zip(self.basis_elements,grouped_tag):
			for tag in tags:
				basis_dict[tag] = elements
		""""""
		qubit_dict = {}
		qubit_index = 0
		for index,atom in enumerate(model):
			if len(basis_dict[atom.tag]) == 1:
				qubits = []
			elif len(basis_dict[atom.tag]) == 2:
				qubits = [qubit_index]
				qubit_index += 1
			elif len(basis_dict[atom.tag]) == 0:
				raise ValueError(f"The tag '{atom.tag}' is not corresponding to any basis.")
			else:
				qubits = [qubit_index+addition for addition in range(len(basis_dict[atom.tag]))]
				qubit_index += len(basis_dict[atom.tag])
			qubit_dict[index]={"variable_num":len(basis_dict[atom.tag]),"basis":basis_dict[atom.tag],"tag":atom.tag,"qubits":qubits}
		return qubit_dict
		## bulkはindexをtagに反映していないかも知れないので、広げた後の対称性を明らかにするために,tagを設定するコードを作る。金属間化合物で検証
		## -> self.tag2grouped_basisがそもそもない slabはbulk2slabでtagつけているし

	def make_QUBO(self,qubit_dict=None,model=None,ECI_file_name=None): #Cythonizeをする必要がある
		## 3次以上で減らせそうなら変数を減らす
		if qubit_dict is None:
			qubit_dict = self.make_qubit(model)
		dimension = (len(self.max_cluster_size)+1) # 4次も3次も次数下げで同じbit数になることを確認済み
		if dimension >= 5:
			# 実はここを消せばもっと高次にも対応しているが、一旦動作確認をplot_QUBOとかでしていないので開放していない
			raise NotImprementedError("This version of SurfQit does not support HUBO for 5 or higher atoms cluster.")
		basis_set         = {} ## tag毎に構築する
		for qubit_info in qubit_dict.values():
			if qubit_info["tag"] not in basis_set.keys():
				basis             = self.get_basis_function(qubit_info["basis"])
				basis_set_unsort  = basis.get_basis_functions()
				basis_set_origin  = []
				for partial_basis_set in basis_set_unsort:
					partial_basis = {}
					for qubit_basis in qubit_info["basis"]:
						partial_basis[qubit_basis] = partial_basis_set[qubit_basis]
					basis_set_origin.append(partial_basis)
				basis_array       = np.array([[phi_i_j for phi_i_j in phi_i.values()] for phi_i in basis_set_origin])
				basis_set[qubit_info["tag"]] = basis_array
		if ECI_file_name is not None:
			cf_names, ECIs    = self.get_ECI(ECI_file_name)
		else:
			cf_names, ECIs    = self.get_ECI()
		cf_dict = {name: 1.0 for name in cf_names}
		if model is not None:
			atoms = model.copy()
		else:
			atoms = self.QUBO_model.copy()
		if self.model_type == "bulk":
			print(f"This bulk system has {len(np.unique(atoms.get_tags()))} types of inequivalent sites, ",*np.unique(atoms.get_tags()))
			print("If this result is different from your intention, please set tag for each atom of the bulk.")
		elif atoms.get_tags().sum() == 0.0:
			raise ValueError("The input model does not have specific tag corresponding to symmetrical equivalent sites.")
		atoms.set_atomic_numbers(atoms.get_tags())
		cluster_index_dict = cluster2index(self.config,self.settings.cluster_mng.prim,self.settings.cluster_list.clusters,atoms)

		if self.annealing_machine == "FA": ## DA, Jijも作る?
			max_qubit_index = 0
			for value in qubit_dict.values():
				max_qubit_index = max(max_qubit_index,max(value["qubits"]))
			try:
				QUBO = np.zeros(tuple([max_qubit_index+1]*dimension)) ## ここからcythonize可能
			except:
				QUBO = np.zeros(tuple([max_qubit_index+1]*dimension),dtype=np.float32) ## ここからcythonize可能, float64だとメモリが足りない場合用。float16に落とすバージョンも作ってもいいかもしれないけど精度...
		constant = 0

		if self.annealing_machine == "FA": ## DA, Jijも作る? Cythonize
			#print("1330",cf_names) # debugで重要。configで[1,1,3.5] などとmax cluster diaを決めたときに、見たくない次数が入っていないことと見たい次数が軽めのcf数で入っていることを確認できる。
			for ECI, cf_name in zip(ECIs,cf_names):
				cluster_name = None
				cluster_indices = None
				for cluster_name_, cluster_indices_ in zip(cluster_index_dict.keys(),cluster_index_dict.values()):
					if cluster_name_ in cf_name:
						cluster_name = cluster_name_
						cluster_indices = cluster_indices_
				if cluster_name == "c0":
					constant += ECI
				elif "c1_" in cluster_name or "c1" == cluster_name: # c11_みたいなものをはじいている。
					basis = basis_set[int(cluster_name.split("_")[-1])][int(cf_name.split("_")[-1])]
					if int(qubit_dict[cluster_indices[0][0]]["tag"]) != int(cluster_name.split("_")[-1]):
						raise ValueError("The tag of the model is invalid.")
					# 対称移動にかかわらず代入する値は同じ(tagは動かない)->先に作れるし、cluster_indices[0]の結果を援用できる
					sub_QUBO     = ECI*(basis_set[int(qubit_dict[cluster_indices[0][0]]["tag"])][int(cf_name.split("_")[-1]),0]-basis_set[int(qubit_dict[cluster_indices[0][0]]["tag"])][int(cf_name.split("_")[-1]),1])/len(cluster_indices)
					sub_constant = ECI*basis_set[int(qubit_dict[cluster_indices[0][0]]["tag"])][int(cf_name.split("_")[-1]),1]/len(cluster_indices)
					sub_onehot   = ECI*basis_set[int(cluster_name.split("_")[-1])][int(cf_name.split("_")[-1]),:]/len(cluster_indices)
					for cluster_index in cluster_indices:
						bit = qubit_dict[cluster_index[0]]['qubits'] # cluster_indexのlenが1なことは確定している。
						if len(bit) == 1: ## 基底が2つ
							QUBO[*tuple(bit*dimension)] += sub_QUBO
							constant  += sub_constant
						else: ## 今はone hot
							for index, bit_num in enumerate(bit):
								QUBO[*tuple([bit_num]*dimension)] += sub_onehot[index]
				else: ##bulkはtagで取り出せない可能性がある。そもそもcf nameの名づけ方が違うし。チェックが必要。
					if self.model_type == "bulk":
						raise NotImplementedError("SurfQit does not support bulk now.")
					basis = None
					for unit_basis in np.meshgrid(*[basis_set[int(tag)][int(basis_index)] for tag,basis_index in zip(cluster_name.split("__")[-1].split("_"),list(cf_name.split("_")[-1]))]):
						if basis is None:
							basis = np.ones_like(unit_basis)
						basis = basis*unit_basis # np.multiplyは挙動が違うので冗長だがこのように書かざるを得ない, c2で01の元素のペアなら01を見ると値が入っている
					trans = list(range(len(basis.shape)))
					trans[0] = 1
					trans[1] = 0
					basis = basis.transpose(*trans)
					dummy_sub   = basis*ECI/len(cluster_indices)
					dummy_bit   = []
					dummy_count = -1
					for index in cluster_indices[0]:
						dummy_list = []
						for qubit in qubit_dict[index]['qubits']:
							dummy_count += 1
							dummy_list.append(str(dummy_count))
						if len(qubit_dict[index]['qubits']) == 1:
							dummy_list.append(str(dummy_count)+"_dummy") # bitとは異なり、0番目が存在(bit=1)に対応。1番目がbit=0に対応。
						dummy_bit.append(dummy_list)
					sub_bit      = [qubit for index in cluster_indices[0] for qubit in qubit_dict[index]['qubits']]
					sub_QUBO     = np.zeros(tuple([len(sub_bit)]*dimension))
					sub_constant = 0
					for dummy_index in itertools.product(*[np.arange(_) for _ in dummy_sub.shape]): # dummy_subは各原子が選ばれた時の代入値が、sub_QUBOは各bitに対応する代入値が、入っている。
						dummy_count = 0
						for count,index in enumerate(dummy_index):
							if "_dummy" in dummy_bit[count][index]: # dummyが入っているとすればlenが1のものにappendで追加されているのでindexは必ず1
								dummy_count += 1
						if dummy_count ==len(dummy_index): # 全部dummy, [0,...]がECIを、他が0を表すようにQUBOを作る
							#sub_QUBO     = np.zeros(tuple([len(sub_bit)]*dimension)) # 検証用
							#sub_constant = 0 # 検証用
							sub_constant += dummy_sub[*dummy_index] #ここ以下は1が入ってくるときにconstantが打ち消されるように(このECIが表現されるように)QUBOを構築
							pair_atom_indices = []
							for count,index in enumerate(dummy_index):
								pair_atom_indices.append(int(dummy_bit[count][index].split("_")[0]))
							pair_atom_dict = {str(dim): [] for dim in range(1,dimension+1)}
							for pair_atom_index in itertools.product(pair_atom_indices,repeat=dimension):
								pair_atom_dict[str(len(set(pair_atom_index)))].append(pair_atom_index)
							"""安定動作版, debugをいろいろしても動くならこれは消す
							pair_atom_dict = {"1":[],"2":[],"3":[],"4":[]} # 今のところ4次までしかサポートしていないのでbitの種類は1~4しかあり得ない
							for pair_atom_index in itertools.product(pair_atom_indices,repeat=dimension):
								if len(set(pair_atom_index)) == 1:
									pair_atom_dict["1"].append(pair_atom_index)
								elif len(set(pair_atom_index)) == 2:
									pair_atom_dict["2"].append(pair_atom_index)
								elif len(set(pair_atom_index)) == 3:
									pair_atom_dict["3"].append(pair_atom_index)
								elif len(set(pair_atom_index)) == 4:
									pair_atom_dict["4"].append(pair_atom_index)
								else:
									raise NotImplementedError("This version of SurfQit is not support more than 5th order binary optimization problem.")
							"""
							coefficients = {}
							for dim in range(1,dimension+1):
								if len(pair_atom_dict[str(dim)])!=0: # 2元系だとpair_atom_dict[3]以降が何も無い(0or1なので)
									pair_atom_index    = pair_atom_dict[str(dim)][0]
									coefficient        = np.zeros(dimension,dtype=np.int16)
									coefficient[dim-1] = 1                                                                               # m^nのmの係数は1に決まっている
									div = dim**len(pair_atom_index)                                                                      # m^nを実際に計算
									for key,value in coefficients.items():                                                               # 下の次数からmと係数ベクトルを取ってくる
										coefficient -= (comb(dim, key, exact=True)*value).astype(np.int16)                               # dimCmに係数をかけて引く(背反を取るため)
										div -= np.sum(((np.arange(dimension)+1)**len(pair_atom_index))*comb(dim, key, exact=True)*value) # 背反を取った後の値を実際に計算, 係数にはそれぞれの底をn乗したものをかける
									coefficients[dim] = coefficient
									for pair_atom_index  in pair_atom_dict[str(dim)]:
										sub_QUBO[*pair_atom_index] += ((-1)**dim)*dummy_sub[*dummy_index]/div
							"""安定動作版, debugをいろいろしても動くならこれは消す""""""
							for pair_atom_index in pair_atom_dict["1"]:
								sub_QUBO[*pair_atom_index] -= dummy_sub[*dummy_index]
							for pair_atom_index in pair_atom_dict["2"]:
								sub_QUBO[*pair_atom_index] += dummy_sub[*dummy_index]/(2**len(pair_atom_index)-2)
							for pair_atom_index in pair_atom_dict["3"]:
								sub_QUBO[*pair_atom_index] -= dummy_sub[*dummy_index]/(3**len(pair_atom_index)-3*(2**len(pair_atom_index))+3)
								""""""
								3^n(n次元)の入力について、
								1種類の係数を取る方法は3C1=3
								2種類の係数を取る方法は3C2=3通りの2種類の決め方があり、各次元にどっちを入れるかなので2^n, ただし、全部の次元がそろうと1種類になってしまうので最初と最後を除いて3*(2^n-2)
								3種類の係数を取る方法は背反なので3^n-3-(3*(2^n-2)) 検算: n=3の時 3! = 6
								さて、dummy_sub[*dummy_index]をaとおくと(各係数において、bitの組合せを考慮した時に)-aを作りたいが
								1種類の係数は3*-a
								2種類の係数は3*(2^n-2)*a/(2^n-2)
								なので、3種類の係数は-a/(3^n-3-(3*(2^n-2)))=-a/(3^n-3*(2^n)+3)
								""""""
							for pair_atom_index in pair_atom_dict["4"]:
								""""""
								4^n(n次元)の入力について、3^n(n次元)の入力のアナロジーで
								1種類の係数を取る方法は4C1=4
								2種類の係数を取る方法は4C2=6通りの2種類の決め方があり、各次元にどっちを入れるかなので2^n, ただし、全部の次元がそろうと1種類になってしまうので最初と最後を除いて6*(2^n-2)
								3種類の係数を取る方法は4C3=4通りの3種類の決め方があり、ダブらせる係数を各次元にどれを入れるかなので3^n, ただし、全部の次元が2つ以下にそろう(3C2*2^n)のはダメ, 1つにそろうときは(1,2の1と1,3の1のように)ダブルカウントされているので3^n-3*2^n+3 検算: n=3の時27個中(1,1,n)が7個、なので(2,2,n)や(3,3,n)とは被りがないことを考えると6が出るはず
								4種類の係数を取る方法は背反なので4^n-4-6*(2^n-2)-4*(3^n-3*2^n+3)=4^n-4-6*2^n+12-4*3^n+12*2^n-12=4^n-4*3^n+6*2^n-4 検算: n=4の時 4! = 24
								さて、dummy_sub[*dummy_index]をaとおくと-aを作りたいが
								1種類の係数は4 * -a = -4a
								2種類の係数は6*(2^n-2) * +a/(2^n-2)=+6a
								3種類の係数は4*(3^n-3*2^n+3) * -a/(3^n-3*(2^n)+3) = -4a
								なので、4種類の係数は+aを作りたくて,+a/(4^n-4*3^n+6*2^n-4)
								""""""
								sub_QUBO[*pair_atom_index] += dummy_sub[*dummy_index]/(4**len(pair_atom_index)-4*(3**len(pair_atom_index))+6*(2**len(pair_atom_index))-4)
								""""""一般化のためのメモ
								アナロジーを利かせると5^nでは
								1種類の係数を取る方法は5C1通り
								2種類の係数を取る方法は5C2通りの2種類の決め方があり、各次元にどっちを入れるかなので2^n, ただし、全部の次元がそろうと1種類になってしまうので最初と最後を除いて10*(2^n-2)
								3種類の係数を取る方法は5C3通りの3種類の決め方があり、ダブらせる係数を各次元にどれを入れるかなので3^n, ただし、全部の次元が2つ以下にそろう(3C2*2^n)のはダメ, 1つにそろうときは(1,2の1と1,3の1のように)ダブルカウントされているので3^n-3*2^n+3 検算: n=3の時27個中(1,1,n)が7個、なので(2,2,n)や(3,3,n)とは被りがないことを考えると6が出るはず
								4種類の係数を取る方法は5C4通りの4種類の決め方があり、ダブらせる係数を各次元にどれを入れるかなので4^n, ただし、全部の次元が3つ以下にそろう(4C3*3^n)のはダメ, 2つ以下にそろうときは(1,2,3の1,2と1,2,4の1,2のように)ダブルカウントされているので+3C2*2^n*2(112か122か), 1つにそろうときは(1,2,3の1, 1,2,4の1, 1,3,4のように)トリプルカウントされているので-2*4C3; 4^n-4*3^n+6*2^n-4
								5種類の係数を取る方法は背反なので5^n-5-10*(2^n-2)-5*(3^n-3*2^n+3)-5*(4^n-4*3^n+6*2^n-4)=5^n-5-10*2^n+20-5*3^n+15*2^n-15-5*4^n+20*3^n-30*2^n+20=5^n-5*4^n+15*3^n-25*2^n+20
								なんか法則性がありそうだが、背反を取り続けるのがよさそう
								"""
						else:
							#sub_QUBO     = np.zeros(tuple([len(sub_bit)]*dimension)) # 検証用
							#sub_constant = 0 # 検証用
							#print(sub_QUBO,sub_constant,dummy_sub[*dummy_index],dummy_bit,dummy_index)
							true_bits  = []
							false_bits = []
							for count,index in enumerate(dummy_index):
								if "_dummy" in dummy_bit[count][index]:
									false_bits.append(int(dummy_bit[count][index].split("_")[0]))
								else:
									true_bits.append(int(dummy_bit[count][index]))
							sub_QUBO[*true_bits+[true_bits[0]]*(dimension-len(true_bits))] += dummy_sub[*dummy_index] # 以下は影響を打ち消し
							for false_length in range(len(false_bits),0,-1):
								for pair_atom_index in itertools.combinations(false_bits,false_length):
									sub_QUBO[*true_bits+[true_bits[0]]*(dimension-false_length-len(true_bits))+list(pair_atom_index)] += ((-1)**false_length)*dummy_sub[*dummy_index]
							"""安定動作版, debugをいろいろしても動くならこれは消す
							if len(false_bits) == 4:
								raise NotImplementedError("This version of SurfQit does not support the 5 =< order QUBO.")
							if len(false_bits) == 3: # 4次までだと全てが1bitで表す場合の0になることはない(上で全て0の場合の処理はしている)。最大が3つ。
								# -1倍でいいのかは不明
								for pair_atom_index in itertools.combinations(false_bits,3):
									print("1456",pair_atom_index)
								sub_QUBO[*true_bits+false_bits] += -dummy_sub[*dummy_index]
							if len(false_bits) >= 2: 
								for pair_atom_index in itertools.combinations(false_bits,2):
									sub_QUBO[*true_bits+[true_bits[0]]*(dimension-2-len(true_bits))+list(pair_atom_index)] += dummy_sub[*dummy_index]
							if len(false_bits) >= 1: 
								for pair_atom_index in itertools.combinations(false_bits,1):
									sub_QUBO[*true_bits+[true_bits[0]]*(dimension-1-len(true_bits))+list(pair_atom_index)] += -dummy_sub[*dummy_index]
							"""
					for cluster_index in cluster_indices:
						bits = list(itertools.chain.from_iterable([qubit_dict[index]['qubits'] for index in cluster_index]))
						constant += sub_constant
						for bit_pair in itertools.product(bits,repeat=len(sub_QUBO.shape)):
							sub_pair = [bits.index(bit) for bit in bit_pair]
							QUBO[*bit_pair] += sub_QUBO[*sub_pair]
		np.save(self.storage_path+"QUBO.npy",QUBO)
		np.save(self.storage_path+"constant.npy",constant)
		return QUBO,constant

import clease.plot_post_process as pp

def plot_fit(evaluator):
	return pp.plot_fit(evaluator)

import clease.plot_post_process as pp
def plot_ECIs(evaluator):
	return pp.plot_eci(evaluator)

import ase.io
import pandas as pd
import matplotlib.pyplot as plt

def plot_QUBO(data_paths,qubit_dict,QUBO,constant,csv_path): # 今はbinaryのみ
	dimension = len(QUBO.shape)
	df = pd.read_csv(csv_path)
	x = []
	y = []
	fig,ax = plt.subplots()
	for path in data_paths:
		bit = np.zeros(QUBO.shape[0],dtype=np.int16)
		model = ase.io.read(path)
		for atom,index in zip(model,qubit_dict.keys()):
			if len(qubit_dict[index]["qubits"])==1: # 単純に[qubit_dict[index]["basis"].index(atom.symbol)]を代入するわけではない。QUBO作成では0側をTrue, 1側をFalseとしているため。ただし、現象論。具体的にどこで変わっているかが今のところよくわからない。
				bit[qubit_dict[index]["qubits"][0]] = abs(1-[qubit_dict[index]["basis"].index(atom.symbol)][0])
			else:
				bit[qubit_dict[index]["qubits"][qubit_dict[index]["basis"].index(atom.symbol)]] = 1
		y_ = np.dot(bit,QUBO)
		for _ in range(dimension-1):
			y_ = np.dot(y_,bit)
		y.append(y_+constant)
		x.append(df[df["id"]==int(path.split("id")[-1].split(".")[0])]["formation_energy"].values[0]/df[df["id"]==int(path.split("id")[-1].split(".")[0])][df.columns[3:]].values.sum())
	ax.plot(x,y,"bo", mfc="none")
	ax.plot([*ax.get_xlim()],[*ax.get_xlim()],c="r")
	ax.plot(x,y,"bo", mfc="none")
	ax.set_title(f"Fit using {len(x)} data points.")
	ax.set_ylabel(r"E$_{QUBO}$ (eV/atom)")
	ax.set_xlabel(r"E$_{DFT}$ (eV/atom)")
	rmse = np.sqrt(np.average((np.array(y)-np.array(x))**2))*1000
	if np.size(x) and np.size(y) != 0:
		e_range = max(np.append(x, y)) - min(np.append(x, y))
		rmin = min(np.append(x, y)) - 0.05 * e_range
		rmax = max(np.append(x, y)) + 0.05 * e_range
	else:
		rmin = -10
		rmax = 10
	ax.axis([rmin, rmax, rmin, rmax])
	ax.text(
		0.95,
		0.01,
		"\n" f"RMSE = {rmse:.3f} meV/atom",
		verticalalignment="bottom",
		horizontalalignment="right",
		transform=ax.transAxes,
		fontsize=12,
	)
	return fig

import numpy as np
from tqdm import tqdm
from amplify import Solver
from amplify.client import FixstarsClient
from amplify import BinarySymbolGenerator, BinaryQuadraticModel
from amplify.constraint import equal_to

class Annealing_Machine_Converter():
	def __init__(self,config,QUBO,constant,qubit_dict):
		self.config          = config
		self.QUBO            = QUBO
		self.constant        = constant
		self.qubit_dict      = qubit_dict
		self.default_penalty = None

	def set_default_penalty(self,reference,method="cost_max"):
		if method == "cost_max":
			self.default_penalty = np.max(np.abs(reference.evaluator.e_dft))
		elif method == "QUBO_max":
			self.default_penalty = np.max(np.abs(reference))*reference.shape[0]
		else: # 源さんに教えていただいた他の方法も実装予定, 目的関数を無くして制約関数だけにして解いてみると、重み係数の良し悪しを排除してどこが課題なのかを抽出することができると思われる。 制約問題として解いた最適解を目的関数のオブジェクトに入れてprint(f.evaluate(result.best.values)を実行し、その結果を上手くスケーリングして合わせると良いと思われる。→より良い方法として、duplicateを用いてその結果の平均値や最大値を使って設定する方法もある。
			raise NotImplementedError(f"This version of SurfQit does not support the this method {method}")

	def bit2energy(self,QUBO,constant,bits):
		energy = np.dot(bits,QUBO)
		for _ in range(len(QUBO.shape)-1):
			energy = np.dot(energy,bits)
		if constant is not None:
			energy += constant
		return energy

	def endcode_to_annealing_machine(self):
		pass

	def setting_constraints(self):
		pass

	def run(self,substitute_number=None,penalty_multipliers=None):
		pass

	def decode_to_universal_style(self, results):
		pass

import numpy as np
from tqdm import tqdm
from amplify import solve, BinarySymbolGenerator, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to,one_hot


class Fixstars_Amplify(Annealing_Machine_Converter):
	def __init__(self,config,QUBO,constant,qubit_dict):
		super().__init__(config,QUBO,constant,qubit_dict)
		self.client = FixstarsClient()
		self.client.parameters.timeout = self.config("annealing","time_out",int) # ms
		self.client.token              = self.config("annealing","token",str)
		self.num_solves                = self.config("annealing","trial",int)
		self.constraints               = []

	def encode_to_annealing_machine(self):
		binary_generator = BinarySymbolGenerator()
		self.bits = binary_generator.array(self.QUBO.shape[0])
		"""過去のバージョン。特にもんだないなら消す
		if len(self.QUBO.shape)==2:
			QUBO = np.dot(np.dot(self.bits,self.QUBO),self.bits.T)
		elif len(self.QUBO.shape)==3:
			QUBO = np.dot(np.dot(self.bits,np.dot(self.bits,self.QUBO)),self.bits)
		elif len(self.QUBO.shape)==4:
			QUBO = np.dot(np.dot(self.bits,np.dot(self.bits,np.dot(self.bits,self.QUBO))),self.bits)
		else:
			raise NotImplementedError("Encoding more than 5D QUBO was not impremented.")
			# for文で処理可能なので将来やる
		"""
		QUBO = self.bit2energy(self.QUBO,None,self.bits)
		return QUBO

	def setting_constraints(self,substitute_dict):
		# TODO: Apply to multi constraints だめなconstraint
		site_dict = {}
		element_dict = {}
		for site,qubit_info in self.qubit_dict.items():
			site_dict[site] = qubit_info["qubits"].copy()
			for element_index,element in enumerate(qubit_info["basis"]):
				if element not in element_dict.keys() and len(qubit_info["basis"]) ==2 and element == qubit_info["basis"][1]:
					element_dict[element] = [[],qubit_info["qubits"].copy()]
				elif element not in element_dict.keys() and len(qubit_info["basis"]) ==2 and element == qubit_info["basis"][0]:
					element_dict[element] = [qubit_info["qubits"].copy(),[]]
				elif element not in element_dict.keys() and len(qubit_info["basis"]) > 2:
					element_dict[element] = [[qubit_info["qubits"].copy()[element_index]],[]]
				elif len(qubit_info["basis"]) ==2 and element == qubit_info["basis"][1]:
					for bit in qubit_info["qubits"].copy():
						element_dict[element][1].append(bit) # bitも1
				elif len(qubit_info["basis"]) ==2 and element == qubit_info["basis"][0]:
					for bit in qubit_info["qubits"].copy():
						element_dict[element][0].append(bit) # bitも0
				elif len(qubit_info["basis"]) > 2:
					element_dict[element][0].append(qubit_info["qubits"].copy()[element_index]) # bitも0
		site_constraint = []
		element_constraint = []
		for site_list in site_dict.values():
			if len(site_list)>1:
				constraint  = one_hot(sum([self.bits[site_index] for site_index in site_list]))
				site_constraint.append(constraint)
		for element, number in substitute_dict.items():
			constraint  = equal_to(sum([1-self.bits[element_index] for element_index in sorted(element_dict[element][1])]+\
											[self.bits[element_index] for element_index in sorted(element_dict[element][0])]),number)
			element_constraint.append(constraint)
		self.constraints = [site_constraint,element_constraint]
		return self.constraints

	def run(self,penalty_multipliers=None,substitute_dict=None):
		QUBO = self.encode_to_annealing_machine()
		if penalty_multipliers is None:
			problem = QUBO
		elif penalty_multipliers is not None and substitute_dict is not None:
			constraints = self.setting_constraints(substitute_dict)
			problem = QUBO
			for multiplier, constraint in zip(penalty_multipliers,constraints):
				if multiplier == 0:
					continue
				if multiplier is None:
					multiplier = self.default_penalty
				if isinstance(constraint,list):
					for same_condition_constraint in constraint:
						problem = problem + multiplier*same_condition_constraint
				else:
					problem = problem + multiplier*constraint
		else:
			raise ValueError("Please set both substitute_number and penalty_multipliers.")
		result = solve(problem,self.client,num_solves=3)
		return result

	def decode_to_universal_style(self,results):
		result_dict = {}
		energy = []
		energy_with_constraints = []
		bits   = []
		for result in results:
			solution = self.bits.decode(result.values)
			energy.append(self.bit2energy(self.QUBO,self.constant,solution))
			energy_with_constraints.append(result.objective)
			bits.append(solution)
		if len(bits)>0:
			result_dict["best_energy"] = np.min(np.array(energy))
			result_dict["best_bits"] = bits[np.argmin(np.array(energy))]
			result_dict["energy"] = energy
			result_dict["energy_with_constraints"] = energy_with_constraints
			result_dict["bits"] = bits
		return result_dict

def best_bits2model(result_dict,ClusterExpansionExcutor):
	model = ClusterExpansionExcutor.QUBO_model.copy()
	for bit_info, atom in zip(ClusterExpansionExcutor.make_qubit().values(),model):
		if len(bit_info['qubits']) > 1:
			atom.symbol = bit_info['basis'][np.argmin((result_dict["best_bits"][bit_info['qubits']]-1)**2)]
		else:
			atom.symbol = bit_info['basis'][((1-result_dict["best_bits"][atom.index])**2).astype(np.int16)]
	return model

	# model.set_atomic_numbers(model.get_tags())
	# cluster_index_dict = cluster2index(config,settings.cluster_mng.prim,cluster_list.clusters,model)
	# bulk2slabの##を確認, 一旦動いたら。
	# bulkでもself.tag2grouped_basis作る?
	# loggerつける?
	# 最終的に11元素以上に対応するためにはcf_nameをcX__tag_tag_tag_basisbasisからcX__tag_tag_tag__basis_basisにする必要はある。ただ今ではない
	"""
	sublattice_list.sort()
	print("474",name_index_relation)
	if "_".join(sublattice_list) not in name_index_relation.keys():
	"""