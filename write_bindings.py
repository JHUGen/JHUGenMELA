enum_dict = {}


with open("/mnt/sda1/RandomRantings/JHU/Andrei_Research/DocumentingMELA/raw_couplings.txt") as f:
    text = f.read()
    for enumcontents in text.split("enum")[1:]:
        enum_size_temp = []
        
        enumcontents = enumcontents.strip()
        # assert enumcontents.startswith("{")
        enum_name = enumcontents.split("{")[0]
        print(f'py::enum_<pymela::{enum_name}>(m, "{enum_name}")')
        enumcontents = enumcontents.split("{")[1].split("}")[0]
        for enumitem in enumcontents.split(","):
            enumitem = enumitem.split("=")[0].strip()
            
            enum_size_temp.append(enumitem)
            
            if "SIZE" in enumitem:
                print(f'\t.value("{enumitem}", pymela::{enumitem});\n')
                for i in enum_size_temp:
                    enum_dict[i] = enumitem
            else:
                print(f'\t.value("{enumitem}", pymela::{enumitem})')
                
with open("/mnt/sda1/RandomRantings/JHU/Andrei_Research/DocumentingMELA/raw_names.txt") as f:
    for line in f:
        var = line.strip().split()[1]
        spinzero = 'nSupportedHiggses' in var
        isInt = 'CLambda' in var
        size_term = ""
        check_size = var.replace(']', '').replace(';', '').split("[")[1:]
        for possible_size in check_size:
            if "SIZE" in possible_size:
                size_term = possible_size
        if check_size:
            var = var[:var.find("[")]
#             string = f"""
# .def("{var}", [](py::object &obj){{
#     Mela &D = obj.cast<Mela&>();
#     return py::array_t<double>(std::vector<int>{{{", ".join(check_size)}}}, (const double*) &D.{var}, obj);
# }})
# """  
            if spinzero:
                if isInt:
                    string = f"MAKE_COUPLING_ARR_SPIN_ZERO({var},{size_term},int)"
                else:
                    string = f"MAKE_COUPLING_ARR_SPIN_ZERO({var},{size_term},double)"
            else:
                string = f"MAKE_COUPLING_ARR_SPIN_ONETWO({var},{size_term})"
            print(string)
        # else:
            # print(var)

with open("/mnt/sda1/RandomRantings/JHU/Andrei_Research/DocumentingMELA/raw_mela_pointers.txt") as f:
    for line in f:
        line = line.strip()
        if "SelfDCoupling" not in line and "SelfDParameter" not in line:
            continue
        elif "SelfDCoupling" in line:
            list_len = 2
        else:
            list_len = 1
        
        
        mela_name, array_name_and_index = line.split("=")
        mela_name = mela_name.strip()

        array_name_and_index = array_name_and_index.replace("SelfDCoupling(", '').replace("SelfDParameter(", '').replace('"', '').replace(")", '').replace("ROOT.pymela.", '')
        array_name_and_index = array_name_and_index.strip().split(',')
        # print(array_name_and_index)
        dimension = len(array_name_and_index)
        
        array_name, index = array_name_and_index[0], array_name_and_index[1:]
        
        if list_len == 2:
            indexstr = ""
            index = [f"[{i.strip()}]" for i in index]
            for i in index:
                indexstr += i
            if len(index) == 2:
                string = f"\t\tMAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO({array_name}, {mela_name}, {index[-1].replace('[', '').replace(']', '')}, 0)"
                print(string)
                
                h_term = mela_name.find('h')
                if h_term == -1:
                    second_higgs_name = mela_name + '_h2'
                else:
                    second_higgs_name = mela_name[:h_term + 1] + '2' + mela_name[h_term + 1:]
                string = f"\t\tMAKE_COUPLING_REAL_IMAGINARY_SPIN_ZERO({array_name}, {second_higgs_name}, {index[-1].replace('[', '').replace(']', '')}, 1)\n"
                print(string)
                
            elif len(index) == 1:
                string = f"\t\tMAKE_COUPLING_REAL_IMAGINARY_SPIN_ONETWO({array_name}, {mela_name}, {index[-1].replace('[', '').replace(']', '')})\n"
                print(string)
        
        elif 'clambda' in line.lower():
            if 'cz_' in mela_name or 'cw_' in mela_name:
                string = f"\t\tMAKE_COUPLING_C_LAMBDA({array_name}, {mela_name}, {index[-1].replace('[', '').replace(']', '')}, 0)"
                print(string)
                
                second_higgs_name = mela_name + "_h2"
                string = f"\t\tMAKE_COUPLING_C_LAMBDA({array_name}, {second_higgs_name}, {index[-1].replace('[', '').replace(']', '')}, 1)\n"
                print(string)
            else:
                string = f"\t\tMAKE_COUPLING_LAMBDA({array_name}, {mela_name}, {index[-2].replace('[', '').replace(']', '')}, {index[-1].replace('[', '').replace(']', '')}, 0)"
                print(string)
                
                second_higgs_name = mela_name + "_h2"
                string = f"\t\tMAKE_COUPLING_LAMBDA({array_name}, {second_higgs_name}, {index[-2].replace('[', '').replace(']', '')}, {index[-1].replace('[', '').replace(']', '')}, 1)\n"
                print(string)
                
                
        #     string = f"""
        # .def_property(
        #     "{mela_name.strip()}", 
        #     py::cpp_function(
        #         [](py::object &obj){{
        #             Mela &D = obj.cast<Mela&>();
        #             return py::array_t<double>(std::vector<int>{{{list_len}}}, (const double*) &D.{array_name}{indexstr}, obj);
        #         }}, py::keep_alive<0, 1>()),
        #     py::cpp_function(
        #         [](Mela &D, std::array<double, 2> coupl){{
        #             D.{array_name}{indexstr}[0] = coupl[0];
        #             D.{array_name}{indexstr}[1] = coupl[1];
        #         }}, py::keep_alive<0, 1>())
        # )"""
            # string = 
        # else:
        #     indexstr = ""
        #     index = [i.strip() for i in index]
        #     indexstr = "{" + ",".join(index) + "}"
        #     string = f"""
        # SELFDPARAMETER.def_property(
        #     "{mela_name.strip()}", 
        #     py::cpp_function(
        #         [](py::object &obj){{
        #             Mela& D = obj.cast<&Mela>();
        #             py::array_t array_val = py::array_t<double>(std::vector<int>{indexstr}, (const double*) &D.{array_name}, obj);
        #             return D.array_val.at({indexstr.replace("{", '').replace("}", '')});
        #         }}),
        #     py::cpp_function(
        #         [](py::object &obj, double coupl){{
        #             Mela &D = obj.cast<Mela&>();
        #             py::array_t array_val = py::array_t<double>(std::vector<int>{indexstr}, (const double*) &D.{array_name}, obj);
        #             array_val.mutable_at({indexstr.replace("{", '').replace("}", '')}) = coupl;
        #         }}, py::keep_alive<0, 1>())
        # )"""