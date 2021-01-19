/******************/
/*--Initializers--*/
/******************/

function Init(){
	d3.tsv(dataLoc+"data/coorProjection.tsv", function(coord){
		coord=coord.map(function(d){
			return {"sample": d.Name, "x" : +d.x, "y": +d.y};
		});
		coordSave=coord;
		xScale = d3.scaleLinear()
		.domain([d3.min(coord, function(d) { return d.x; })*xCoef , d3.max(coord, function(d) { return d.x; })*xCoef])
		.range([padding, w-padding]);
		yScale = d3.scaleLinear()
		.domain([d3.max(coord, function(d) { return d.y; }), d3.min(coord, function(d) { return d.y; })])
		.range([padding, h-padding]);
		cellPoint=svg.selectAll("circle")
		.data(coord)
		.enter()
		.append("circle")
		.attr("cx", function(d) {
			return xScale(d.x);
		})
		.attr("cy", function(d) {
			return yScale(d.y);
		})
		.attr("fill","white")
		.attr("stroke","black")
		.attr("stroke-width",strokeScale(radiusRange.property("value")))
		.attr("r",radiusScale(radiusRange.property("value")))
		.attr("opacity",0.8)
		.on("click",function(d){
			d3.event.stopPropagation();
			selectedSelAnnot=d[selectableAnnotation];
			SelAnnotSelect.select("option[value='"+selectedSelAnnot+"']").property("selected",true);
			selectSelAnnot();
		});
		initColor();
	});
}

function initColor(){
	d3.json(dataLoc+"data/colorScales.json", function(data) {
		predifinedOrdinalColors=data;
		getSampleAnnot();
		initBranchgeneClust();
	});
}

function getSampleAnnot(){
	d3.tsv(dataLoc+"data/sampleAnnot.tsv", function(annot){
		sampleAnnot=annot;
		curData=copyData(sampleAnnot);
		SelAnnotList=curData.map(getSelAnnot).unique();
		SelAnnotSelect
			.on("change",function(){
				selectedSelAnnot=SelAnnotSelect.property("value");
				selectSelAnnot();
			})
			.selectAll("option")
			.data(SelAnnotList)
			.enter()
			.append("option")
			.property("value",function(d){
				return d;
			})
			.text(function(d){
				return d;
			})
		SelAnnotSelect.append("option")
			.property("hidden",true)
			.property("selected",true)
			.property("id","defaultSelAnnotOption")
			.text("Select "+selectableAnnotation);
		cellPoint
			.data(curData)
			.attr("SelAnnot",function(d){
				return d[selectableAnnotation];
			});
		selectSampleCor.selectAll("option")
			.data(curData)
			.enter()
			.append("option")
			.property("value",function(d){
				return d.Name;
			})
			.text(function(d){
				return d.Name;
			})
		initToolTip();
	});
}


function initBranchgeneClust(){
	d3.json(dataLoc+"data/GenesClustList.json",function(geneClustList){
		selectgeneClust.selectAll("option")
			.data(geneClustList)
			.enter()
			.append("option")
			.property("value",function(d){
				return d;
			})
			.text(function(d){
				return d;
			})
	});
	d3.tsv(dataLoc+"data/GenesClust.tsv",function(geneClustData){
		featureData=geneClustData;
		featureData.map(function(d){d.Membership=+d.Membership});
		featureData=featureData.sort(function(x, y){
			return (y.Membership - x.Membership);
		})
	});
}

function initToolTip(){
	var tooltipColorScale=new colorScale(tooltipBackgroundAnnot);
	
	cellPoint
	.data(curData)
	.on("mouseover", function(d) {		
		tipTable.html("");
		tooltip.transition()		
			.duration(transLen/5)		
			.style("opacity", .9);		
		
		Object.entries(d).forEach(function(entry){tipTable.append("tr").html("<td><b>"+entry[0]+":</b></td><td>"+entry[1]+"</td>")})	
		tooltip.style("background",tooltipColorScale.scale(d[tooltipBackgroundAnnot]))
			.style("color",whiteOrBlack(tooltipColorScale.scale(d[tooltipBackgroundAnnot])))
			.style("left", (d3.event.pageX) + "px")		
			.style("top", (d3.event.pageY) + "px");
	})					
	.on("mouseout", function(d) {		
		tooltip.transition()		
			.duration(transLen/5)		
			.style("opacity", 0);
	});
	getGeneList();
}

function getGeneList(){
	d3.json(dataLoc+"data/geneList.json", function(data) {
		geneList=data;
		$( "#searchGeneText" ).autocomplete({
			maxResults: 25,
			source: function(request, response) {
				var results = $.ui.autocomplete.filter(geneList, request.term);
				response(results.slice(0, this.options.maxResults));
			}
		})
		colorize();
	});
}

/******************/
/* --- Colors --- */
/******************/

var colorScale=function(name,type,domain,range){
	var ordinal=true;
	var scale;
	if(name in predifinedOrdinalColors.range){
		domain=predifinedOrdinalColors.domain[name];
		range=predifinedOrdinalColors.range[name];
		scale=d3.scaleOrdinal()
		.domain(domain)
		.range(range);
	}else if(name in sampleAnnot[0]){
		sampleAnnot.map(function(d){d[name]=+d[name]})
		domain=[d3.min(sampleAnnot, function(d) { return d[name]; }) , d3.max(sampleAnnot, function(d) { return d[name]; })]
		range=["white","black"];
		ordinal=false;
		scale=d3.scaleLinear()
		.domain(domain)
		.range(range);
	}
	else{
		if(typeof type === 'undefined' || typeof domain === 'undefined' || typeof range === 'undefined'){
			throw "Missing parameters";
		}
		if(type=="ordinal"){
			scale=d3.scaleOrdinal().domain(domain).range(range);
		}else if(type=="linear"){
			var ordinal=false;
			scale=d3.scaleLinear().domain(domain).range(range);
		}else{
			throw "scale type must be 'ordinal' or 'linear'";
		}
	}
	
	JsonScale=[];
	for(i=0;i<domain.length;i++) JsonScale.push({"domain":domain[i],"range":range[i]});
	
	this.JsonScale=JsonScale;
	this.scale=scale;
	this.ordinal=ordinal;
	this.name=name;
	this.range=range;
	this.domain=domain;
	this.length=domain.length;
	this.maxSlider=100;
	if(ordinal) destroySlider();
	else{
		this.min=domain[0];
		this.max=domain[domain.length-1];
		this.scale100=d3.scaleLinear().domain([0,100]).range([this.min,this.max]);
		this.metrics=calcMetrics(curData.map(function(d){ return d[name]  }));
		if(this.name=="expression" && (! logScaleCheck.property("checked"))){
			this.maxSlider=Math.max(d3.quantile(curData.map(function(d){ return d[name]  }).sort((a, b) => a - b),0.75)/this.max*100,2); //minimum 2%
		}
	}
	this.setLinearDomain=function(newDomain){
		var oldDomain = this.domain;
		if(this.length==2){
			this.domain=newDomain;
		}else{
			this.domain=[newDomain[0],this.domain[1],newDomain[1]];
		}
		
		this.scale=d3.scaleLinear().domain(this.domain).range(range);
		legendSVG
			.selectAll("text")
			.data(this.domain)
			.text(function(d) {
				return(round3(d))
			});		
		updateColor();
	}
}
	


function colorize(){
	destroyViolin();
	destroyTabList();
	exprOnly.style("display","none");
	curData=copyData(sampleAnnot);
	switch(MainAccordion.accordion( "option", "active" )){
		case 0:
			colorizeSampleAnnot();
		break;
		case 1:
			colorizeGene();
		break;
		case 2:
			colorizegeneClust();
		break;
		case 3:
			colorizeSampleCor();
		break;
	}	
}

function colorizeSampleAnnot(){
	var name=selColorAnnot.property("value");
	legendName = selColorAnnot.select("option[value='"+name+"']").text();
	curColorScale=new colorScale(name);
	changeColor();
}

function colorizeGene(){
	var geneName = textGene.property("value");
	legendName = geneName + " normalized counts"
	exprOnly.style("display","block");
	if(geneList.includes(geneName)){
		textGene.classed("is-invalid",false);
		d3.tsv(dataLoc+"data/gene/"+geneName+".tsv",function(gene){
			var data = gene.map(function(d){ return +d.expression });
			if(logScaleCheck.property("checked")) data=data.map(function(d){return Math.log2(d+1)});
			for(i=0; i<gene.length; i++){
				curData[i]["expression"]= data[i];
			}
			curColorScale=new colorScale("expression","linear",[0 , d3.max(data)],["white","red"] );
			createGeneTabList(dataLoc+"data/geneCor/"+geneName+".tsv",["Gene","Pearson correlation"]);		
			changeColor();
			
			var module="";
			featureData.map(function(x){if(x.Name==geneName) module = x.Module})
			if(module != ""){
				geneClustlink.html("<a><h2>geneClust module: "+module+"</h2></a>")
					.on("click",function(d){
					selectgeneClust.property("value",module);
					MainAccordion.accordion("option","active", 2);
			})
			}else{
				geneClustlink.html("");
			}
			
		});
	}else{
		if(geneName==""){
			textGene.classed("is-invalid",false);
		}else{
			textGene.classed("is-invalid",true);
		}
		neutralColorize();
	}
}


function colorizeSampleCor(){
	var branch = selectSampleCor.property("value");
	legendName = "Correlation with "+branch
	d3.tsv(dataLoc+"data/corSamples.tsv",function(SampleCor){
		for(i=0; i<SampleCor.length; i++) curData[i]["Correlation"]=+SampleCor[i][branch];
		var range = [d3.min(curData, function(d) { return d.Correlation; }) ,d3.median(curData, function(d) { return d.Correlation; }) ,d3.max(curData, function(d) { return d.Correlation; })];
		curColorScale=new colorScale("Correlation","linear",range,["#00008B","white","orange"] );
		changeColor();
	});
}

function colorizegeneClust(){
	var module = selectgeneClust.property("value");
	legendName = module + " cluster eigengene"
	d3.tsv(dataLoc+"data/GenesClustEigen.tsv",function(moduleEigen){
		for(i=0; i<moduleEigen.length; i++) curData[i]["Eigen"]=+moduleEigen[i][module];
		var range = [d3.min(curData, function(d) { return d.Eigen; }) ,d3.median(curData, function(d) { return d.Eigen; }) ,d3.max(curData, function(d) { return d.Eigen; })];
		curColorScale=new colorScale("Eigen","linear",range,["blue","white","red"] );
		creategeneClustTablist(module);
		changeColor();
	});
}

function createGeneTabList(file){
	d3.tsv(file,function(gene){
		rightBottomContainer
			.style("opacity",1)
			.style("display","block");
		
		headRightBottomContainer.text("Best/worst Pearson correlations");
		tablegene.html("");
		tablegene.append("thead").append("tr").html("<th>Gene</th><th>Correlation</th>");
		tablegene
			.append("tbody")
			.selectAll("tr")
			.data(gene)
			.enter()
			.append("tr")
			.on("click",function(d){
				textGene.property("value",d.Gene);
				colorize();
			})
			.html(function(d){
				return "<td>"+d.Gene+"</a></td><td>"+round3(d.pearsonCor)+"</td>"
			})			
	});
}

function creategeneClustTablist(module){
	rightBottomContainer
		.style("opacity",1)
		.style("display","block");
			
	headRightBottomContainer.text("Genes of this cluster");
	moduleData=featureData.filter(function(x){return x.Module==module});
	moduleData
	
	tablegene.html("");
	tablegene.append("thead").append("tr").html("<th>Gene</th><th>Module Membership</th>");

	tablegene
		.append("tbody")
		.selectAll("tr")
		.data(moduleData)
		.enter()
		.append("tr")
		.on("click",function(d){
			textGene.property("value",d.Name);
			MainAccordion.accordion("option","active", 1);
		})
		.html(function(d){
			return "<td>"+d.Name+"</a></td><td>"+round3(+d.Membership)+"</td>"
		})			
}


function destroyTabList(){
	headRightBottomContainer.text("");
	rightBottomContainer
		.style("opacity",0)
		.style("display","none");
}

function changeColor(){
	traceLegend();
	if(!curColorScale.ordinal) plotViolin();
	updateColor();
}

function changeViolin(){
	destroyViolin();
	plotViolin();
}

function updateColor(){
	cellPoint
	   .data(curData)	
	   .attr("fill", function(d) {
			return curColorScale.scale(d[curColorScale.name]);
	   });	
}

function neutralColorize(){
	curColorScale=null;
	destroySlider();
	destroyViolin();
	d3.selectAll(".exprOnly").style("display","none")
	legendName= " ";
	cellPoint
	   .data(sampleAnnot)
	   .attr("fill", "white");	
	legendTitle
			.transition()		
			.duration(transLen/2)
			.style("opacity", 0)
			.on("end",function(){
				legendTitle.text(" ");
			})
	legendSVG
			.transition().duration(transLen/2)
			.attr("height",0)
			.on("end",function(){
				legendSVG.html("")
			});
}

/**********************/
/*--DOM manipulation--*/
/**********************/

function traceLegend(){
	legendTitle
			.transition()		
			.duration(transLen/2)
			.style("opacity", 0)
			.on("end",function(){
				legendTitle.html(legendName);
			})
			.transition()
			.duration(transLen/2)
			.style("opacity", 1);
			
	legendSVG
		.transition().duration(transLen/2)
		.attr("height",0)
		.on("end",function(){
			legendSVG.html("")
			.selectAll("rect")
				.data(curColorScale.JsonScale)
				.enter()
				.append("rect")
				.attr("width",squareDim)
				.attr("height",squareDim)
				.attr("stroke","black")
				.attr("stroke-width",1)
				.attr("y",function(d,i){
					return 1+i*(squareDim+padding)
				})
				.attr("x",1)
				.attr("fill",function(d){
					return d.range;
				});
			
			legendSVG
				.selectAll("text")
				.data(curColorScale.JsonScale)
				.enter()
				.append("text")
				.text(function(d){
					return curColorScale.ordinal?d.domain:round3(d.domain);
				})
				.attr("font-size",fontSize)
				.attr("x",(squareDim+padding))
				.attr("y",function(d,i){
					return (squareDim+fontSize)/2.3+i*(squareDim+padding)
				});
			legendSVG
				.transition().duration(transLen/2)
				.on("start",function(){
					if(!curColorScale.ordinal){
							$("#slider").slider({
								range: true,
								values: [0,curColorScale.maxSlider]
							})
							.slider({
							  slide: function( event, ui ) {
								curColorScale.setLinearDomain(ui.values.map(function(d){
									return curColorScale.scale100(d)
								}));
							  }
						});
						sliderContainer.style("visibility","visible")
							.transition().duration(transLen/3)
							.style("opacity", 1)
						curColorScale.setLinearDomain($("#slider").slider("values").map(function(d){
							return curColorScale.scale100(d)
						}));
					}
				})
				.attr("height",(squareDim+padding)*curColorScale.length+1);
		});
}


function selectSelAnnot(){
	nodeSelAnnot=svg.selectAll("circle[SelAnnot='"+selectedSelAnnot+"']");
	nodeNotSelAnnot=svg.selectAll("circle").filter(function(d){return d[selectableAnnotation]!=selectedSelAnnot});
	nodeSelAnnot
		.attr("stroke-width", strokeScale(radiusRange.property("value")*1.75))
		.attr("r",radiusScale(radiusRange.property("value")*1.75))
		.attr("opacity", 0.9)
	nodeNotSelAnnot
		.attr("stroke-width", strokeScale(radiusRange.property("value")))
		.attr("opacity", 0.1)
		.attr("r",radiusScale(radiusRange.property("value")))
}

function unselectSelAnnot(){
	selectedSelAnnot=null;
	nodeSelAnnot=null;
	nodeNotSelAnnot=null;
	d3.select("#defaultSelAnnotOption").property("selected",true)
	cellPoint
		.attr("stroke","black")
		.attr("stroke-width",strokeScale(radiusRange.property("value")))
		.attr("opacity",0.8)
		.attr("r",radiusScale(radiusRange.property("value")));
}

function changeRadiusWTanim(){
	if(selectedSelAnnot===null){
		cellPoint
		.attr("r",radiusScale(radiusRange.property("value")))
		.attr("stroke-width",strokeScale(radiusRange.property("value")));
	}else{
		nodeSelAnnot
			.attr("r",radiusScale(radiusRange.property("value"))*1.75)
			.attr("stroke-width",strokeScale(radiusRange.property("value"))*1.75);
		nodeNotSelAnnot
			.attr("r",radiusScale(radiusRange.property("value")))
			.attr("stroke-width",strokeScale(radiusRange.property("value")));
	}
}

function plotViolin(){
		colLegend.style("border-bottom","2px solid #D8D8CB");
		quantScaleOnly.style("display","block");
		var margins = {top: 15, right: 50, bottom: 120, left: 75};
		var yViolin = violinSelect.property("value");
		WdivDistro = divDistro.node().getBoundingClientRect().width;
		HdivDistro = WdivDistro*9/16;
		chart = makeDistroChart({
            data:curData,
            xName:yViolin,
            yName:curColorScale.name,
			scale: "linear",
			margin:margins,
            axisLabels: {xAxis: mainAnnot, yAxis: legendName},
            selector:"#chart-distro",
            chartSize:{
				width:WdivDistro, 
				height:HdivDistro
			},
            constrainExtremes:true,
			Xdomain: predifinedOrdinalColors.domain[yViolin],
			colors:predifinedOrdinalColors.range[yViolin]
		});
		chart.renderViolinPlot({
			width: 80,
			resolution: 50,
			clamp: 0
		});
		chart.renderBoxPlot({
			showWhiskers:false,
			showOutliers:false,
			boxWidth:10,
			lineWidth:10,
			colors:['#555']
		});
		d3.selectAll(".chart-wrapper .box-plot .box").attr("fill-opacity",0.4).attr("stroke-width",2);
		d3.selectAll(".chart-wrapper .box-plot line").attr("stroke-width",2);
		d3.selectAll(".chart-wrapper .box-plot circle").attr("fill","white").attr("stroke","black");
		d3.selectAll(".chart-wrapper .box-plot .median").attr("stroke","black");
		d3.selectAll(".chart-wrapper .box-plot circle.median").attr("fill","white");
		d3.selectAll(".chart-wrapper .box-plot .mean").attr("stroke","white").attr("stroke-width",1).attr("stroke-dasharray","2,1");
		d3.selectAll(".chart-wrapper .violin-plot .area").attr("opacity",0.8).attr("shape-rendering","geometricPrecision");
		d3.selectAll(".chart-wrapper .violin-plot .line").attr("stroke-width",0.75).attr("fill","none").attr("shape-rendering","geometricPrecision");
		d3.selectAll(".chart-wrapper text").attr("font-family","Arial").attr("font-size","1em");
		d3.selectAll(".chart-wrapper .axis path, .chart-wrapper .axis line")
			.attr("fill","none")
			.attr("stroke","#888")
			.attr("stroke-width",2)
			.attr("shape-rendering","crispEdges");
		d3.selectAll(".chart-wrapper .y.axis .tick line")
			.attr("stroke","black")
			.attr("opacity",0.6)
			.attr("stroke-dasharray","2,1")
			.attr("stroke-width","1")
			.attr("shape-rendering","crispEdges");
		
		d3.select(".chart-wrapper .backgroundRect").attr("fill","#D8D8CB")
		d3.select("#chart-distro svg").property("export","violinPlot").classed("exportable",true)
		
		var wViolinLegend= violinLegend.node().width.baseVal.value
		var hViolinLegend= wViolinLegend*svgRatio;
		var xScaleViolinLegend = d3.scaleLinear()
			.domain([d3.min(coordSave, function(d) { return d.x; })*xCoef , d3.max(coordSave, function(d) { return d.x; })*xCoef])
			.range([padding, wViolinLegend-padding]);
		var yScaleViolinLegend = d3.scaleLinear()
			.domain([d3.max(coordSave, function(d) { return d.y; }), d3.min(coordSave, function(d) { return d.y; })])
			.range([padding, hViolinLegend-padding]);
		
		var violinLegendPoint=violinLegend.selectAll("circle")
		.data(coordSave)
		.enter()
		.append("circle")
		.attr("cx", function(d) {
			return xScaleViolinLegend(d.x);
		})
		.attr("cy", function(d) {
			return yScaleViolinLegend(d.y);
		})
		.attr("stroke","black")
		.attr("stroke-width",0.2)
		.attr("r",2.5)
		.attr("opacity",0.8)
		
		violinLegendPoint.data(curData).attr("fill",function(d){
			return chart.colorFunct(d[yViolin]);
		});
}
	
	
function destroyViolin(){
	d3.selectAll(".quantScaleOnly").style("display","none")
	divDistro.html("");
	colLegend.style("border-bottom","none");
	violinLegend.html("");
	chart=null;
}

var lineGen = d3.line()
    .x(function(d) { return xScale(d.x); })
    .y(function(d) { return yScale(d.y); });
	
var radiusScale=d3.scaleLinear()
	.domain([1,100])
	.range([0.3,10]);
var strokeScale=d3.scaleLinear()
	.domain([1,100])
	.range([.05,1.5]);
	
function resizeSVG(){
	w = svg.node().getBoundingClientRect().width;
	h=w*svgRatio;
	svg.attr("width",w)
		.attr("height", h)
		.call(zoom);
}

function destroySlider(){
	if(sliderContainer.style("visibility")=="visible"){
		sliderContainer.transition().duration(transLen/3)
			.style("opacity", 0)
			.on("end",function(){
				sliderContainer.style("visibility","hidden");
				$("#slider").slider("destroy");
			});		
	}
}



function saveSvg() {
    d3.selectAll(".exportable").each(function(){
		var element=d3.select(this)
		element.node().setAttribute("xmlns", "http://www.w3.org/2000/svg");
		var svgData = element.node().outerHTML;
		var preface = '<?xml version="1.0" standalone="no"?>\r\n';
		var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
		var svgUrl = URL.createObjectURL(svgBlob);
		var downloadLink = document.createElement("a");
		downloadLink.href = svgUrl;
		downloadLink.download = element.property("export")+".svg";
		document.body.appendChild(downloadLink);
		downloadLink.click();
		document.body.removeChild(downloadLink);
	});
   
}

function downloadGeneTable(){
	var breakChar="\r\n"
	var separator=","
	out=geneTabListContainer.select("thead").selectAll("th").nodes()
		.map(function(x){return x.textContent}).join(separator);
	out+=breakChar;
	geneTabListContainer.select("tbody").selectAll("tr").each(function(d,i){ 
		out+=d3.select(this).selectAll("td").nodes().map(function(y){return y.textContent}).join(separator)
		out+=breakChar;
	});
	
	var textBlob = new Blob([out], {type: "text/csv;charset=utf-8"});
	var textUrl = URL.createObjectURL(textBlob);
	var downloadLink = document.createElement("a");
	downloadLink.href = textUrl;
	downloadLink.download = "PTUI_tableExport.csv";
	document.body.appendChild(downloadLink);
	downloadLink.click();
	document.body.removeChild(downloadLink);
}


/*****************/
/*-- Utilities --*/
/*****************/
var getSelAnnot=function(that){return that[selectableAnnotation]};

function round3(x){
	return Math.round(1000*x)/1000;
}

function copyData(oldData){
	var newData=[];
	oldData.forEach(function(d){
		newData.push(Object.assign({}, d));
	});
	return newData;
}

Array.prototype.unique = function() {
  return this.filter(function (value, index, self) { 
    return self.indexOf(value) === index;
  });
}


function whiteOrBlack(color) {

    // Variables for red, green, blue values
    var r, g, b, hsp;
    
    // Check the format of the color, HEX or RGB?
    if (color.match(/^rgb/)) {

        // If HEX --> store the red, green, blue values in separate variables
        color = color.match(/^rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*(\d+(?:\.\d+)?))?\)$/);
        
        r = color[1];
        g = color[2];
        b = color[3];
    } 
    else {
        
        // If RGB --> Convert it to HEX: http://gist.github.com/983661
        color = +("0x" + color.slice(1).replace( 
        color.length < 5 && /./g, '$&$&'));

        r = color >> 16;
        g = color >> 8 & 255;
        b = color & 255;
    }
    
    // HSP (Highly Sensitive Poo) equation from http://alienryderflex.com/hsp.html
    hsp = Math.sqrt(
    0.299 * (r * r) +
    0.587 * (g * g) +
    0.114 * (b * b)
    );

    // Using the HSP value, determine whether the color is light or dark
	if (hsp>150) {

        return '#000000';
    } 
    else {

        return '#FFFFFF';
    }
}