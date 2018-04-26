function new_struct=sortStructByValue(new_struct, value)
[unused, order] = sort([old_struct(:).value]);
new_struct = old_struct(order);