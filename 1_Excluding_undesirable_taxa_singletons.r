###

# Set colors for plotting
phylum_colors <- c(
  "gray",'black', "#DA5724", "#5F7FC7","#508578", "#CD9BCD", "orange",
  "#5F7FC7","#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#D1A33D", "#8A7C64", "#599861","#5E738F"
)

my_color_collection <- c(
  "#CBD588", "#5F7FC7", "orange", "#AD6F3B", "#673770",
  "#D14285", "#652926", "#C84248", "#8569D5", "#5E738F",
  "#D1A33D", "#8A7C64", "#599861","#616163", "#FFCDB2",
  "#6D9F71", "#242F40",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")

my_color_Class <- c(
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "orange", "#5F7FC7", "#CBD588", "#AD6F3B", "#673770",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',

  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")



my_color_Order <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

my_color_Family <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)
my_color_Family2 <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

# up to 300
my_color_OTU <- c(
  "gray",'black', "#5F7FC7", "orange",  "#AD6F3B",
  "#673770","#D14285", "#652926", "#C84248",  "#8569D5",
  "#5E738F","#D1A33D", "#8A7C64", "#599861","#616163",
  "#FFCDB2", "#242F40", "#6D9F71",  "#CCA43B", "#F92A82",
  "#ED7B84", "#7EB77F", "#DEC4A1", "#E5D1D0", '#0E8482',
  '#C9DAEA', '#337357', '#95C623', '#E55812', '#04471C',
  '#F2D7EE', '#D3BCC0', '#A5668B', '#69306D', 'navy',
  '#1A535C', '#4ECDC4', 'orange', '#FF6B6B', "orchid1",
  'cyan2', '#FFF275', 'springgreen', '#FF3C38', '#A23E48',
  '#000000', '#CF5C36', '#EEE5E9', '#7C7C7C', '#EFC88B',

  '#2E5266', '#6E8898', '#9FB1BC', '#D3D0CB', '#E2C044',
  '#5BC0EB', '#FDE74C', '#9BC53D', '#E55934', '#FA7921',
  "#CD9BCD", "#508578", "#CBD588","#CBD588", "#5F7FC7",
  "orange",   "#AD6F3B", "#673770","#D14285", "#652926",
  "#C84248",  "#8569D5", "#5E738F","#D1A33D", "#8A7C64",
  "#599861","#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", "#DEC4A1",
  "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', '#95C623',
  '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', '#A5668B',
  '#69306D', '#0E103D', '#1A535C', '#4ECDC4', '#F7FFF7',
  '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange")

my_color_otu2 <- c(
  "gray",'black', "#AD6F3B", "#673770","#D14285",
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D",
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40",
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB', 'pink',
  "orange", "#CBD588", "#5F7FC7",
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "orange")


my_color_gen <- c(
  "white",'yellow', "#AD6F3B", "#673770","#D14285",
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D",
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40",
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#CD9BCD', '#6699CC', 'pink',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB',

  "orange", "#CBD588", "#5F7FC7",
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357')


### Bacteria

#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
phy.cp <- subset_taxa(phy, Class == "D_2__Chloroplast") ## just confirming
phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 130 otus of CP
vec.cp

# (2) MT
#phy.mt <- subset_taxa(phy, Family == "D_4__Mitochondria")
phy.mt <- subset_taxa(phy, Order == "D_3__Rickettsiales") #337
vec.mt <- rownames(otu_table(phy.mt))
tax_table(phy.mt)
length(rownames(otu_table(phy.mt))) ## 337 otus of CP

# (3) Unassigned
unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy, Kingdom == "Unassigned")
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 4 otus

### exclude those vectors
phy  #17133 taxa and 1176 samples

sample_names(phy)
sample_variables(phy)

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

phy.clean <- pop_taxa(phy, c(vec.cp, vec.mt, vec.un))

##before clean-up
phy  ## 28901
sum(otu_table(phy)) ##25067648

##after clean-up
phy.clean ##16662
sum(otu_table(phy.clean)) # 12239207

# checking procedure of whether the MT and CP otus are cleaned
taxa_names(subset_taxa(phy, Order == "D_3__Rickettsiales"))  # before 337
taxa_names(subset_taxa(phy, Order=="D_3__Chloroplast")) # before 130

taxa_names(subset_taxa(phy.clean , Order == "D_3__Rickettsiales")) # after 0
taxa_names(subset_taxa(phy.clean , Family=="D_4__Mitochondria")) # after 0
taxa_names(subset_taxa(phy.clean , Order == "D_3__Chloroplast")) # after 0



#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt()
# phy.test4

tax_table(phy.clean)[,colnames(tax_table(phy.clean))] <- gsub(tax_table(phy.clean)[,colnames(tax_table(phy.clean))],pattern="[A-Z]_[0-9]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(phy.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

phy.clean    ## 28330
str(phy.clean)
otu_table(phy.clean)

phy.clean  ## 28330


# ## fix it in phy.clean object!!! pop_taxa does the work
# phy.clean <- pop_taxa(phy.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(phy.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement
phy.clean.otu <- otu_table(phy.clean)
head(phy.clean.otu)
df.clean.otu <- data.frame(phy.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


sample_names(phy.clean)



## Remove reads with over 300 bp
library(seqinr)
bac.seq <- read.fasta(file = "dna_sequences_bac.fasta", as.string = TRUE, seqtype = "DNA")
bac.seq[which(getLength(bac.seq)>300)]
otu_over_300bp <- attr(bac.seq[which(getLength(bac.seq)>300)], "names")


#whole
#2018
phy.clean.ss
phy.clean.ss <- pop_taxa(phy.clean,otu_over_300bp)
phy.clean.ss  ## 16605
sum(otu_table(phy.clean.ss)) # 12197391


## Divide bacteria and archaea
# #2017
# arch.clean.ss.17 <- subset_taxa(phy.clean.ss.17, Kingdom == "Archaea")
# arch.clean.ss.17 <- phyloseq::filter_taxa(arch.clean.ss.17, function(x) sum(x) != 0, TRUE)
#
# bac.clean.ss.17 <- subset_taxa(phy.clean.ss.17, Kingdom != "Archaea")
# bac.clean.ss.17 <- phyloseq::filter_taxa(bac.clean.ss.17, function(x) sum(x) != 0, TRUE)
#
# #2018
# arch.clean.ss.18 <- subset_taxa(phy.clean.ss.18, Kingdom == "Archaea")
# arch.clean.ss.18 <- phyloseq::filter_taxa(arch.clean.ss.18, function(x) sum(x) != 0, TRUE)
#
# bac.clean.ss.18 <- subset_taxa(phy.clean.ss.18, Kingdom != "Archaea")
# bac.clean.ss.18 <- phyloseq::filter_taxa(bac.clean.ss.18, function(x) sum(x) != 0, TRUE)


### whole
arch.clean.ss <- subset_taxa(phy.clean.ss, Kingdom == "Archaea")
arch.clean.ss <- phyloseq::filter_taxa(arch.clean.ss, function(x) sum(x) != 0, TRUE)

bac.clean.ss <- subset_taxa(phy.clean.ss, Kingdom != "Archaea")
bac.clean.ss <- phyloseq::filter_taxa(bac.clean.ss, function(x) sum(x) != 0, TRUE)


###Fungal community
#### get unassigned vectors

# (3) Unassigned
unique(tax_table(fun)[,'Kingdom']) ## "k__Chromista", "k__Plantae", "Unassigned"
tax_table(fun)[,'Kingdom']
fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Chromista","k__Plantae","k__Rhizaria","k__unidentified"))
vec.un <- rownames(otu_table(fun.un))
tax_table(fun.un)
length(rownames(otu_table(fun.un))) ##  191

### exclude those vectors
fun  # 328 samples, 10981 taxa

sample_names(fun)
sample_variables(fun)

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### Clean it up!!!
fun.clean <- pop_taxa(fun, c(vec.un))

#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(fun.two)[,colnames(tax_table(fun.two))] <- gsub(tax_table(fun.two)[,colnames(tax_table(fun.two))],pattern="[A-Z]_[0-9]__",replacement="")
# fun.test4 <- fun.two %>% psmelt()
# fun.test4
tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean    ## 5968 taxa
str(fun.clean)
otu_table(fun.clean)

fun.clean  ## 10865


# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement
fun.clean.otu <- otu_table(fun.clean)
head(fun.clean.otu)
df.clean.otu <- data.frame(fun.clean.otu)
dim(df.clean.otu)
df.clean.otu$total <- apply(df.clean.otu, 1, sum)
head(df.clean.otu)
df.clean.otu <- tibble::rownames_to_column(df.clean.otu, var = 'OTU')


sample_names(fun.clean)

##### get rid of otu of less than 100 reads

fun.seq <- read.fasta(file = "dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
fun.seq[which(getLength(fun.seq)<100)]
#
otu_less_than_100bp <- attr(fun.seq[which(getLength(fun.seq)<100)], "names")
fun.clean.ss
fun.clean.ss <- pop_taxa(fun.clean,otu_less_than_100bp)
fun.clean.ss ## 9437
sum(otu_table(fun.clean.ss)) ## 35113746

## Fungal OTUs which are 'unidentifed' at the phylum level
tax.fun<-tax_table(fun.clean.ss)
tax.fun <-data.frame(tax.fun)
unidentified.phylum<-rownames(tax.fun)[is.na(tax.fun$Phylum)]

unidentified.phylum.seq<-fun.seq[which(attr(fun.seq, "names")%in% unidentified.phylum)]

write.fasta(unidentified.phylum.seq, names(unidentified.phylum.seq), file.out = "unidentified.phylum.seq.fasta")
length(names(unidentified.phylum.seq))

## Remove unidentified OTUs based on BLAST result
unmatched
remove.fungal.otu <- names(unidentified.phylum.seq)[-which(names(unidentified.phylum.seq)%in% unmatched)]

fun.clean.ss2 <- pop_taxa(fun.clean.ss, remove.fungal.otu)



## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id

bac.list
fun.list

otu.list <- rbind(bac.list, fun.list)
dim(otu.list)
# write.xlsx(otu.list,'otu_id_book.xlsx')


OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],arch.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id
