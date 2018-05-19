'''Module to parse input argument to obtain the ontology label field number(s)'''

from .ProcLog import ProcLog
import re

'''
Function name: _convert_to_ontology_label_index
Inputs       : Expression specifying ontology label field numbers
Outputs      : Ontology label field numbers
Description  : Converts multiple values or value ranges into a set of ontology field numbers.
               For example, '2,4,-1--3' is converted to a set (-1, -2, -3, 2, 4)
'''
def _convert_to_ontology_label_index(otl_label_field_num_opt):
    ontology_label_field_indices = set()
    label_field_ranges = otl_label_field_num_opt.split(',')
    for label_field_range in label_field_ranges:
        if re.match(r'^-?\d+$', label_field_range):
            val = int(label_field_range)
            if val > 0:
                val -= 1
            elif val == 0:
                ProcLog.log_exec_error('Ontology label field number cannot be zero')
                break
                
            ontology_label_field_indices.add(val)
        else:
            m = re.match(r'^(-?\d+)-(-?\d+)$', label_field_range)
            if m:
                val1 = int(m.group(1))
                val2 = int(m.group(2))
                if val1 > 0:
                    val1 -= 1
                    
                if val2 > 0:
                    val2 -= 1

                if val1 > val2:
                    ontology_label_field_indices.update(range(val2, val1 + 1))
                else:
                    ontology_label_field_indices.update(range(val1, val2 + 1))
                    
    return ontology_label_field_indices

'''
Function name: parse_ontology_label_field_nums
Inputs       : User's argument specifying ontology label field numbers
Outputs      : Extracted ontology label field numbers
Description  : First replaces the '~' symbol by '-' to represent negative values (as '-' is not
               allowed in input argument), and then converts the input argument to field numbers
'''
def parse_ontology_label_field_nums(otl_label_field_num_opt):
    otl_label_field_num_opt = otl_label_field_num_opt.replace('~', '-')

    if re.match(r'^(-?\d+(--?\d+)?,)*(-?\d+(--?\d+)?)$', otl_label_field_num_opt):
        return _convert_to_ontology_label_index(otl_label_field_num_opt)
    else:
        ProcLog.log_exec_error('Invalid ontology label field specification: {}'.format(otl_label_field_num_opt))
        return None

'''test'''
'''
ProcLog.init_logs()
otl_indices = parse_ontology_label_num_fields('2,4-5')
print(sorted(otl_indices))
otl_indices = parse_ontology_label_num_fields('4-7')
print(sorted(otl_indices))
otl_indices = parse_ontology_label_num_fields('4-7,2,9')
print(sorted(otl_indices))
otl_indices = parse_ontology_label_num_fields('1,4-3,6-7')
print(sorted(otl_indices))
otl_indices = parse_ontology_label_num_fields('-1--3,2,3-5,8')
print(sorted(otl_indices))
'''
