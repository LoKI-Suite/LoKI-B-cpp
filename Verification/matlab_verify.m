input_path = '~/github/luxurious-loki/Verification/input/';
type_path = 'growth/';

for n=100:100:10000
	display([type_path num2str(n) '.in']);
	lokibcl([input_path type_path num2str(n) '.in']);
end

type_path = 'ee/';

for n=1000:1000:10000
	display([type_path num2str(n) '.in']);
	lokibcl([input_path type_path num2str(n) '.in']);
end	
