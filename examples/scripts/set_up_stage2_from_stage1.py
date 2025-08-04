import argparse
import os
import xml.etree.ElementTree as ET


# parent -> children
def get_parent2children_tag() -> dict[str, list[str]]:
    return {
        'rawReadCountsModel': ['adoRate'],
        'seqCovModel': ['allelicSeqCov', 'allelicSeqCovRawVar'],
        'nucReadCountsModel': ['effSeqErrRate', 'shapeCtrl1', 'shapeCtrl2'],
        'siteModel': ['shape']
    }


def get_manipulated_bg_num_tag() -> dict[str, str]:
    return {
        'bgNum': 'bgSitesNum',
        'bg2varRatio': 'bgSites2VarSitesRatio'
    }


# parameter -> (lower, upper)
def get_param_bounds() -> dict[str, tuple[float|None, float|None]]:
    return {
        'shape': (6.0, None)
    }


# tag -> attribute
def get_tree_index() -> dict[str, str]:
    return {
        'tree': 'treeFileName'
    }


def get_name2id(
        parent2children_tag: dict[str, list[str]],
        param_bounds: dict[str, tuple[float|None, float|None]],
        root: ET.Element
) -> dict[str, tuple[str, tuple[float|None, float|None]]]:
    name2id: dict[str, tuple[str, tuple[float|None, float|None]]] = {}
    for parent, children in parent2children_tag.items():
        for parent_tag in root.iter(parent):
            for attr_name, attr_val in parent_tag.attrib.items():
                if attr_name in children and attr_val.strip().startswith('@'):
                    __attr_val = attr_val.strip().lstrip('@')
                    if attr_name in param_bounds.keys():
                        name2id[attr_name] = (__attr_val, *(param_bounds[attr_name]))
                    else:
                        name2id[attr_name] = (__attr_val, None, None)

    for parent_tag in root.iter("state"):
        for child in parent_tag:
            if child.tag == "parameter":
                attr_val = child.attrib['id']
                if attr_val not in name2id.keys() and attr_val in param_bounds.keys():
                    name2id[attr_val] = (attr_val, *(param_bounds[attr_val]))
                else:
                    name2id[attr_val] = (attr_val, None, None)

    return name2id


def get_estimate_type(results: str, key: str) -> None|str:
    if os.path.isfile(results):
        with open(results, 'r') as fh:
            for line in fh:
                if key in line.strip():
                    if 'mean' in line:
                        return 'mean'
                    elif 'median' in line:
                        return 'median'
                    else:
                        return 'mode_gaussian'


def get_estimates(
        estimates: str,
        name2id: dict[str, tuple[str, tuple[float|None, float|None]]],
        estimate_type: str
) -> None|dict[str, float]:
    if os.path.isfile(estimates):
        id2estimates: dict[str, float] = {}
        reach_beginning: bool = False
        reach_ending: bool = False
        with open(estimates, 'r') as fh:
            for line in fh:
                if reach_ending is True:
                    break

                if reach_beginning is True:
                    if not line.strip().startswith('#'):
                        comp: list[str] = line.strip().split('\t')
                        if comp[1].strip() in estimate_type:
                            for _name, _collect in name2id.items():
                                _id, _lower, _upper = _collect
                                if comp[0].strip() == _id:
                                    _val = float(comp[2].strip())
                                    if (_lower is not None and _lower > _val) or (_upper is not None and _upper < _val):
                                        raise ValueError(f'Error! Parameter {_id} ({_val}) is out of bounds: [{_lower}, {_upper})]. Handel it manually.')
                                    id2estimates[comp[0].strip()] = _val
                    else:
                        reach_beginning = False
                        reach_ending = True
                else:
                    if line.startswith('#MCMC samples'):
                        reach_beginning = True
            return id2estimates


def update_template2(
        tree2: ET.ElementTree,
        parent2children_tag: dict[str, list[str]],
        name2id: dict[str, tuple[str, float|None, float|None]],
        id2estimates: dict[str, float],
        tree_index: dict[str, str],
        tree_file_name: str|None,
        manipulated_bg_tag: dict[str, str],
        bg_num: float|None = None,
        bg2var_ratio: float|None = None
) -> None:
    for parent, children in parent2children_tag.items():
        counter: int = 0
        for parent_tag in tree2.getroot().iter(parent):
            counter += 1
            if counter > 1:
                raise ValueError(f'Error! More than 1 tag with name {parent} found, but only one allowed.')
            for child_tag in parent_tag.iter():
                if child_tag.tag in name2id.keys() or (child_tag.tag == 'parameter' and child_tag.get('name') in name2id.keys()):
                    if child_tag.tag in name2id.keys():
                        key = child_tag.tag
                    else:
                        key = child_tag.get('name')
                    child_tag.text = str(id2estimates[name2id[key][0]])

    if tree_file_name is not None:
        for _tag, _attr in tree_index.items():
            for tree_tag in tree2.getroot().iter(_tag):
                if _attr in tree_tag.attrib.keys():
                    tree_tag.attrib[_attr] = os.path.abspath(tree_file_name)

    for parent_tag in tree2.getroot().iter("state"):
        for child in parent_tag:
            if child.tag == "parameter" and child.attrib['id'] in id2estimates.keys():
                child.text = str(id2estimates[child.attrib['id']])
    
    if len(tree2.getroot().findall('data')) != 1:
        raise ValueError(f'Error! Multiple `data` tags found for the target template. Only one allowed.')
    data_tag: ET.Element = tree2.getroot().findall('data')[0]

    if bg_num is not None and bg_num > 0:
        for attr in manipulated_bg_tag.values():
            if attr in data_tag.attrib.keys():
                data_tag.attrib.pop(attr)
        data_tag.set(manipulated_bg_tag['bgNum'], str(bg_num))

    if bg2var_ratio is not None and bg2var_ratio > 0:
        for attr in manipulated_bg_tag.values():
            if attr in data_tag.attrib.keys():
                data_tag.attrib.pop(attr)
        data_tag.set(manipulated_bg_tag['bg2varRatio'], str(bg2var_ratio))


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(
        description='Setup configuration files for stage 2 based on results of stage 1.'
    )

    sys_args_parser.add_argument(
        '--tree',
        type=str,
        required=False,
        help='starting tree for stage 2 with the MCC tree from stage 1'
    )
    sys_args_parser.add_argument(
        '--estimates',
        type=str,
        required=True,
        help='estimated parameter values from stage 1'
    )
    sys_args_parser.add_argument(
        '--results',
        type=str,
        required=True,
        help='a file containing the a list of files of variants called from stage 1'
    )
    sys_args_parser.add_argument(
        '--template1',
        type=str,
        required=True,
        help='template configuration file for stage 1'
    )
    sys_args_parser.add_argument(
        '--template2',
        type=str,
        required=True,
        help='template configuration file for stage 2'
    )
    sys_args_parser.add_argument(
        '--bgNum',
        type=int,
        default=None,
        help='the number of background sites'
    )
    sys_args_parser.add_argument(
        '--bg2varRatio',
        type=float,
        default=None,
        help='the ratio between the number of background sites and the number of variant sites'
    )
    sys_args_parser.add_argument(
        '--out',
        type=str,
        required=True,
        help='modified template configuration file for stage 2'
    )

    sys_args = sys_args_parser.parse_args()
    if sys_args.tree is not None and not os.path.isfile(sys_args.tree):
        raise ValueError(f'Error! {sys_args.tree} does not exist.')
    if not os.path.isfile(sys_args.estimates):
        raise ValueError(f'Error! {sys_args.estimates} does not exist.')
    if not os.path.isfile(sys_args.results):
        raise ValueError(f'Error! {sys_args.results} does not exist.')
    if not os.path.isfile(sys_args.template1):
        raise ValueError(f'Error! {sys_args.template1} does not exist.')
    if not os.path.isfile(sys_args.template2):
        raise ValueError(f'Error! {sys_args.template2} does not exist.')
    if sys_args.bgNum is not None and sys_args.bg2varRatio is not None:
        raise ValueError(f'Error! `bgNum` and `bg2varRatio` are mutually exclusive.')
    if not os.path.isdir(os.path.dirname(sys_args.out)):
        os.makedirs(os.path.dirname(sys_args.out), exist_ok=True)
    return sys_args


def main():
    parent2children_tag = get_parent2children_tag()
    param_bounds = get_param_bounds()
    tree_index = get_tree_index()

    sys_args = parse_system_args()

    tree1 = ET.parse(
        sys_args.template1,
        ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    )

    name2id = get_name2id(parent2children_tag, param_bounds, tree1.getroot())
    estimate_type = get_estimate_type(sys_args.results, '.vcf')
    id2estimates = get_estimates(
        sys_args.estimates,
        name2id,
        estimate_type
    )

    tree2 = ET.parse(
        sys_args.template2,
        ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    )
    update_template2(
        tree2,
        parent2children_tag,
        name2id,
        id2estimates,
        tree_index,
        sys_args.tree,
        get_manipulated_bg_num_tag(),
        sys_args.bgNum,
        sys_args.bg2varRatio
    )
    tree2.write(open(sys_args.out, 'w'), encoding='unicode')


if __name__ == '__main__':
    main()
