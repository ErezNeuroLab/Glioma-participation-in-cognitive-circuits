% create DAN cog data

DAN_heart = add_cog_data(1, 'heart');
DAN_salt_lighthouse = add_cog_data(25, 'salt_lighthouse');
DAN_salt_RT = add_cog_data(27, 'salt_RT');

DAN_table = readtable('DAN_table.txt');
DAN_cog_data = [DAN_table, DAN_heart, DAN_salt_lighthouse, DAN_salt_RT];

writetable(DAN_cog_data);