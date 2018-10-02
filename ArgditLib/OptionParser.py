'''Module to parse input argument to obtain the sequence class label field number(s)'''

from .ProcLog import ProcLog
import re

'''
Function name: _convert_to_seq_class_label_index
Inputs       : Expression specifying sequence class label field numbers
Outputs      : Sequence class label field numbers
Description  : Converts multiple values or value ranges into a set of class label field
               numbers. For example, '2,4,-1--3' is converted to a set (-1, -2, -3, 2, 4)
'''
def _convert_to_seq_class_label_index(class_label_field_num_opt):
    class_label_field_indices = set()
    label_field_ranges = class_label_field_num_opt.split(',')
    for label_field_range in label_field_ranges:
        if re.match(r'^-?\d+$', label_field_range):
            val = int(label_field_range)
            if val > 0:
                val -= 1
            elif val == 0:
                ProcLog.log_exec_error('Sequence class label field number cannot be zero')
                break
                
            class_label_field_indices.add(val)
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
                    class_label_field_indices.update(range(val2, val1 + 1))
                else:
                    class_label_field_indices.update(range(val1, val2 + 1))
                    
    return class_label_field_indices

'''
Function name: parse_seq_class_label_field_nums
Inputs       : User's argument specifying sequence class label field numbers
Outputs      : Extracted sequence class label field numbers
Description  : First replaces the '~' symbol by '-' to represent negative values (as '-' is not
               allowed in input argument), and then converts the input argument to field numbers
'''
def parse_seq_class_label_field_nums(class_label_field_num_opt):
    class_label_field_num_opt = class_label_field_num_opt.replace('~', '-')

    if re.match(r'^(-?\d+(--?\d+)?,)*(-?\d+(--?\d+)?)$', class_label_field_num_opt):
        return _convert_to_seq_class_label_index(class_label_field_num_opt)
    else:
        err_msg = 'Invalid sequence class label field specification: {}'.format(class_label_field_num_opt)
        ProcLog.log_exec_error(err_msg)
        return None

'''test'''
'''
ProcLog.init_logs()
class_label_indices = parse_seq_class_label_field_nums('2,4-5')
print(sorted(class_label_indices))
class_label_indices = parse_seq_class_label_field_nums('4-7')
print(sorted(class_label_indices))
class_label_indices = parse_seq_class_label_field_nums('4-7,2,9')
print(sorted(class_label_indices))
class_label_indices = parse_seq_class_label_field_nums('1,4-3,6-7')
print(sorted(class_label_indices))
class_label_indices = parse_seq_class_label_field_nums('-1--3,2,3-5,8')
print(sorted(class_label_indices))
'''
