    margin = {top: 50, right: 20, bottom: 50, left: 20};
	Colors = {};
	Colors.names=getColors();

 
	//override window.open function
  window.open =function(open){
    return function(url,name,features){
      var testPop=open("","","width=100,height=100");
      if (testPop===null){
        alert("please disable your popup blocker");
      }
      testPop.close();
      name=name||"blank";
      return open.call(window,url,name,features);
    };
  }(window.open);
    	
  checkPopBlocker=function(){
    var testPop=window.open("","","width=100,height=100");
    if (testPop===null){
      alert("please disable your popup blocker");
    }
    testPop.close();
  }


	function onload(){
		
		s1=[];  //selected group1 
		s2=[];  //selected group2 
		currentSelected=[]; //current selected 
		
		//------------------------------------------------------------
		RL=parseJSON(data);
		nodes=buildTree(RL[0],RL[1]);
		
		colorList=getColors();
		colorList=Object.values(colorList);
		
		if (typeof drdata!=="undefined"){
			vizcells=JSON.parse(drdata);
		}
		
		root=nodes[0];
		cells=RL[2];
		edges=RL[1];
		GL=RL[3].map(function(d){return d.toUpperCase();});
		TFs=RL[4];
		ETFs=RL[5];
		drawTree(root);
		drawCells(cells);
		drawEdgeTexts(TFs);
		updateNodes(nodes);
		settextcolor();

		// load TF-DNA data
		geneRegMap=JSON.parse(tfdna);
		regTargetMap=buildTargetMap(geneRegMap); //reg -> gene targets

		//------------------------------------------------------------
		
		// create dropdown menu
		
		createDropDown("Edge",'#edgeGroup1','g1edgeid',edges.map(function(d){return 'P'+d.id+":"+'N'+nodes[d.source_node].id+'->'+'N'+nodes[d.end_node].id;}),function(){return gedgeonchange(1);})
		createDropDown("Edge",'#edgeGroup2','g2edgeid',edges.map(function(d){return 'P'+d.id+":"+'N'+nodes[d.source_node].id+'->'+'N'+nodes[d.end_node].id;}),function(){return gedgeonchange(2)})
		
		// select nodes
		var ds=new DragSelect({
				selectables: document.getElementsByClassName("cells"),
				multiSelectKeys: ['shiftKey'],
				callback: function(elements){
					for (var e of elements){
						if (window.event.shiftKey){
							d3.selectAll(elements)
							.style("fill","#aaa")
							currentSelected.push(e);
						}
					}
				  
				}
			});		
	}
	
	///////////////////////////////////////////////////////////////////
	
	
	
	//expression in selected
	function exSelect(){
		var inputsgene=document.getElementById("scellGeneName").value;
		if (s1.length==0 && s2.length==0){
			alert("please select cells!")
		}else{
            cellids=cells.map(function(d){return d.id;});
            s1ids=s1.map(function(d){return cellids.indexOf(d.__data__.id);});
            s2ids=s2.map(function(d){return cellids.indexOf(d.__data__.id);});
			$.get("/ht-bin/colorCellsC.py", {"cellGeneName":inputsgene},function(data,status){
				var exdata=JSON.parse(data);
                if (s1ids.length>0){
                    var ex1=exdata.filter(function(d,i){
                        if (s1ids.indexOf(i)!=-1){
                            return true;
                        }return false;
                    });

                    var ex1sum=ex1.reduce(function(a,b){return a+b;});
                    ex1sum=ex1sum/ex1.length;
                }
                if (s2ids.length>0){
                    var ex2=exdata.filter(function(d,i){
                        if (s2ids.indexOf(i)!=-1){
                            return true;
                        }return false;
                    });

                   
                    var ex2sum=ex2.reduce(function(a,b){return a+b;});
                    ex2sum=ex2sum/ex2.length;

                }
				
                var tsum=exdata.reduce(function(a,b){return a+b;});
                if (ex1sum==undefined){
                    var ex1sum=(tsum-ex2sum)/(exdata.length-ex2.length);
                    barplot(["unselected","group2"],[ex1sum,ex2sum],"expression of "+inputsgene+" in selected cells");
                }

                else if (ex2sum==undefined){
                    var ex2sum=(tsum-ex1sum)/(exdata.length-ex1.length);
                    barplot(["group1","unselected"],[ex1sum,ex2sum],"expression of "+inputsgene+" in selected cells");
                }

                else if (ex1sum!=undefined && ex2sum!=undefined){
                    barplot(["group1","group2"],[ex1sum,ex2sum],"expression of "+inputsgene+" in selected cells");
                }
                
			});
		}
			
	}
	
	//showSpring
	
	function ShowSpring(){
		if (s1.length==0 && s2.length==0){
			alert("please select cells!")
		}else{
			
			var SCors=[];
			var SGroups=[];
			var SLabels=[];
			var CellIDs=cells.map(function(d){return d.id;});
			var s1IDs=s1.map(function(d){return d.__data__.id;});
			var s2IDs=s2.map(function(d){return d.__data__.id;});
			var scolor1=document.getElementById("scolor1").value;
			var scolor2=document.getElementById("scolor2").value;
			
			
			var ucolor="#D3D3D3"
			for (var icor of CellCor){
				var cid=icor[0];
				if (s1IDs.indexOf(cid)!=-1){
					SCors.push([parseFloat(icor[1]),-1*parseFloat(icor[2])])
					SGroups.push("Group1");
					SLabels.push(scolor1);
					
				}
				else if(s2IDs.indexOf(cid)!=-1){
					SCors.push([parseFloat(icor[1]),-1*parseFloat(icor[2])])
					SGroups.push("Group2");
					SLabels.push(scolor2);
				}else if (CellIDs.indexOf(cid)!=-1){
					SCors.push([parseFloat(icor[1]),-1*parseFloat(icor[2])]);
					SGroups.push("Unselected");
					SLabels.push(ucolor);
				}
				
			}
			//console.log(SLabels);
			scatterplot(SCors,SGroups,SLabels);
			window.open("https://kleintools.hms.harvard.edu/tools/springViewer_1_6_dev.html?cgi-bin/client_datasets/nacho_springplot/allMerged","","width=800,height=800");
		}
		
	}



	function compareSelect(){

		if (s1.length==0 && s2.length==0){
			alert("please select cells!")
		} else{
			var cellsID=cells.map(function(d){return d.id;});
			var scells1=s1.map(function(d){return cellsID.indexOf(d.__data__.id);});
			var scells1=scells1.join(",");
			var scells2=s2.map(function(d){return cellsID.indexOf(d.__data__.id);});
			var scells2=scells2.join(",");

			var nw=window.open("","","width=600,height=600");
			nw.document.write("<html><body><p id='runtag'>running...</p></body></html>");
			var tdiv = d3.select(nw.document.body)
			.append("div");

			$.post("/ht-bin/deSelect.py",{"cells1":scells1,"cells2":scells2},function(data,status){
				var deList=JSON.parse(data);
				//console.log(deList);
				createTable(tdiv,"detable",[["Gene","Group2","Group1","log2 fold change","Transcription factors regulating the gene"]].concat(deList));
				nw.document.getElementById("runtag").innerHTML="done!";
				//

				tdiv
				.append("p")
				.text("Table: Differentially expressed genes (UP-regulated) between Group1 and  Group2");

				var species="MOUSE"
				//add GO analysis
				tdiv
				.append("button")
				.text("functional analysis")
				.on("click",function(){
					panthergoInput(deList.map(function(d){return d[0];}), species);
				});
				//

			});
		}
	}
	
	//clear selected
	function clearSelect(){
		s1=[];
		s2=[];
		currentSelected=[];
		var cellabels=[];
		for (var ci of cells){
			if (cellabels.indexOf(ci.type)==-1){
				cellabels.push(ci.type);
			}

		}
		
		var usedcolors=[];
		for (var li in cellabels){
			var lii=cellabels[li];
			var clii=Colors.names[Object.keys(Colors.names)[li]];
			usedcolors.push([lii,clii]);
		}
		d3.selectAll("circle")
		.style("fill",function(d){
			var ci=cellabels.indexOf(d.type);
			var ci=ci%Object.keys(Colors.names).length;
			var ckey=Object.keys(Colors.names)[ci];
			var cl=Colors.names[ckey];
			return cl;
			})
		.style("opacity",0.7);
		document.getElementById("sgroup1").innerHTML=s1.length;
		document.getElementById("sgroup2").innerHTML=s2.length;
		
	}
	

	//set selected
	function setSelect(groupid){
		var gcolor=document.getElementById("scolor"+groupid).value;
		d3.selectAll(currentSelected)
		.style("fill",gcolor);
		document.getElementById("sgroup"+groupid).innerHTML=currentSelected.length;
		if (groupid==1){
			s1=currentSelected;
		}else{
			s2=currentSelected;
		}
		currentSelected=[];
	}

	//

	function ptest(){
			document.getElementById("cellStatus").innerHTML="coloring...";
			var inputgene=document.getElementById("cellGeneName").value;
			var excut=parseInt(document.getElementById("cellGeneSlider").innerHTML);
			$.get("/ht-bin/colorCellsC.py", {"cellGeneName":inputgene},function(data,status){
				var exdata=JSON.parse(data);
				var cmax=Math.max.apply(Math,exdata);
				var cmin=Math.min.apply(Math,exdata);

				var myColor = d3.scale.linear().domain([cmin,cmax]).range(["lightgray", "red"])

				d3.selectAll(".cells")
				.style("fill",function(d,i){
					return myColor(exdata[i]);
					})
                .attr("r",function(d,i){
                    if (exdata[i]>excut){
                        return 5;
                    } else {
                    	if (exdata[i]> 0) {
                    		return 3;
						} else {
                    		return 1;
						}
					}
                })
                d3.selectAll(".cells")
                .style("opacity",function(d,i){
                    if (exdata[i]>0){
                        return 1;
                    } return 0.2;
                })
				document.getElementById("cellStatus").innerHTML="coloring done";

			});

	}


	//this function is used to pre-compute all needed data for nodes (e.g. target fold change).
	function updateNodes(nodes){
		for (var node of nodes){
			nodeUpdate(node);
		}	
	}
	
	//calculate fold change for all gene for given node (parent->node)
	//calculate edge ending at chosen node 

	function nodeUpdate(node){
		var pnode=node.parent;
		var chosenEdge;
		var fc={};
		for (var edge of edges){
			if (nodes[edge.to]==node){
				chosenEdge=JSON.parse(JSON.stringify(edge));
			}
		}
		if (chosenEdge!=null){
			for (var i in GL){
				fcg=node.E[i]-pnode.E[i];
				gi=GL[i];
				fc[gi]=fcg;
			}
		}
		node["fc"]=fc;
		node["endEdge"]=chosenEdge;
	}

	function resetconfig(){
		zoom(50);
		document.getElementById('zoomsliderbar').value=50;
		document.getElementById("bgcolor").value="#333333";
		setbgcolor();
		resetPath();
		location.reload();
	}
	
	
	//set bgcolor
	function setbgcolor(){
		var color=document.getElementById("bgcolor").value;
		d3.select("svg").style("background", color);
	}

	//set nodecolor
	function setnodecolor(){
		var nodecolor=document.getElementById("nodecolor").value;

		d3.selectAll("g.node")
		.selectAll("circle")
		.attr("fill",nodecolor)
		.attr("stroke",nodecolor)
		.attr("stoke-width","3px")
	}
	
	//set textcolor
	function settextcolor(){
		var textcolor=document.getElementById("textcolor").value;

		d3.selectAll("text")
		.style("fill",textcolor);
		
		
	}
	
	function setlinecolor(){
		var linecolor=document.getElementById("linecolor").value;
		d3.selectAll("line")
		.style("stroke",linecolor)
	}
	
	
	//handle explore tf

	function exploretf(){
		resetPath();
		var tfinput=d3.select("#tfName").property("value").toUpperCase();
		var paths=d3.select("svg").selectAll("path")[0];
		var selectedpaths=[];
		var tnodeList=[];
		var etfList=[];
		for (var edge of edges){
			var etf=edge.etf.map(function(d){return d[1];});
			if(etf.indexOf(tfinput)!=-1){
				tnodeList.push(nodes[edge.to]);
				etfList.push(etf.indexOf(tfinput));
			}
		}
		
		var tpathList=[];
		for (var path of paths){
			var pt=path.__data__.target;
			if (tnodeList.indexOf(pt)!=-1){
				index_pt=tnodeList.indexOf(pt);
				tpathList.push([path,etfList[index_pt]]);
			}
		}
			
		
		for (var path of tpathList){
			var colorpath=hsl_col_perc(path[1]/60,0,120);
			d3.select(path[0]).attr("stroke",colorpath);
		}	
		
		
		//path text
		var textpaths=d3.select("svg").selectAll(".linktext")[0];
		var tpathIndex=tpathList.map(function(d){return paths.indexOf(d[0]);});
		var ttextpath=tpathIndex.map(function(d){return textpaths[d];});
		for (i in ttextpath){
			var tfi=etfList[i];
			var tpath=ttextpath[i];
			d3.select(tpath).text("p-value rank:"+tfi)
			.attr("fill","#fff");
		}
		if (GL.indexOf(tfinput)!=-1){		
			plottf(tfinput);
		}
		//alert("implements explore tf:(1) highlight tf re-regualting path; (2): plot tf expression (3) plot tf target expression, (4) more");
	}
	
	
	function scatterplot(xData,xLabels,xTypeLabels){
			var dlabels={};
			var clabels={};
			for (var ix in xData){
				var di=xData[ix];
				var li=xLabels[ix];
				if (!(li in dlabels)){
					dlabels[li]=[{"x": di[0],"y": di[1]}];
					clabels[li]=xTypeLabels[ix];
				}else{
					dlabels[li].push({"x": di[0],"y": di[1]});
				}
			}
				
			var CombinedDataSet=[];
			var colorCounter=0;
			var dlabels_keys=Object.keys(dlabels).sort();
			for (var di in dlabels_keys){
				var xlabel=dlabels_keys[di];
				var xList=dlabels[xlabel];
				var clabel=clabels[xlabel];
				//CombinedDataSet.push({"label": xlabel, "data": xList, "backgroundColor": Object.values(Colors.names)[colorCounter]});
				CombinedDataSet.push({"label": xlabel, "data": xList, "backgroundColor": clabel, "pointRadius": 1});
				colorCounter+=1;
			}
			
			var plotdata={"datasets":CombinedDataSet};
			
			var newW3 = open('','_blank','height=800,width=800')
			newW3.document.write('<head><title>scatter plot</title> </head><body></body>');
			d3.select(newW3.document.body).append("canvas")
						.attr("id","scatterplot")
						.attr("width","800px")
						.attr("height","800px");
						
			var ctx=newW3.document.getElementById("scatterplot").getContext('2d');
			ctx.globalAlpha=0.4;
			
			var scatterChart = new Chart(ctx, {
				type: 'scatter',
				data: plotdata,
				options: {
				responsive : false,
				scales: {
					xAxes: [{
						type: 'linear',
						position: 'bottom',
						ticks:{
							min: -2500,
							max: 2500
						}
					}],
					yAxes:[{
						ticks:{
							reverse:true,
							min: -2500,
							max: 2500
							}
					}]
				},
				tooltips: {
					callbacks:{
						label: function(tooltipItem,data){
							var pos=[tooltipItem.xLabel,tooltipItem.yLabel];
							for (var ix in  xData){
								if (xData[ix].toString()==pos.toString()){
									return xTypeLabels[ix];
								}							
							}
						}
					}
				}
			}
				
			});	
		}
	
	function barplot(X,Y,Z,P){
			//X: x-axis values
			//Y: y-axis values
			//Z: legend
			
			var newW3 = open('','_blank','height=600,width=800,left='+P);
			newW3.document.write('<head><title></title> </head><body></body>');
			d3.select(newW3.document.body).append("canvas")
						.attr("id","lineplot")
						.attr("width","700px")
						.attr("height","500px");
						
			var ctx=newW3.document.getElementById("lineplot").getContext('2d');
			
			var lineChart = new Chart(ctx, {
				type: 'bar',
				data:{
					labels: X,
					datasets:  [{
						data: Y,
						backgroundColor: "blue",
						label: ""
					}]
				},
				options: {
					responsive : false,
					title: {
						display: true,
						text: Z
					},
                    scales:{
                        yAxes:[{
                            ticks:{
                                beginAtZero:true
                            }
                        }]
                    }
				
				}
			});
		}

	//plot tf
	function plottf(tfinput){
		var LX=nodes.map(function(d){return 'N'+d.id;});
		var LY=nodes.map(function(d){return (d.expression[GL.indexOf(tfinput)]).toFixed(3);});
		barplot(LX,LY,"Expression of "+tfinput);
	
	}

	function plottf_edges(tfinput){
		var LX=edges.map(function(d){return 'P'+d.id;});
		var LY=edges.map(function(d){return (d.expression[GL.indexOf(tfinput)]).toFixed(3);});
		barplot(LX,LY,"Expression of "+tfinput);
	}

	//handle explorede
	function explorede(){
		resetPath();
		//alert("this");
		var deinput=d3.select("#deNameNode").property("value").toUpperCase();
		
		/*
		var paths=d3.select("svg").selectAll("path")[0];
		var selectedpaths=[];
		var tnodeList=[];
		var deList=[];
		for (var edge of edges){
			var de=edge.de.map(function(d){return d.toUpperCase();});
			if(de.indexOf(deinput)!=-1){
				tnodeList.push(nodes[edge.to]);
				deList.push(de.indexOf(deinput)/de.length);
			}
		}
		
		var tpathList=[];
		for (var path of paths){
			var pt=path.__data__.target;
			if (tnodeList.indexOf(pt)!=-1){
				index_pt=tnodeList.indexOf(pt);
				tpathList.push([path,deList[index_pt]]);
			}
		}
			
		
		
		for (var path of tpathList){
			var colorpath=hsl_col_perc(path[1],0,120);
			d3.select(path[0]).attr("stroke",colorpath);
		}	
		
		//path text
		var textpaths=d3.select("svg").selectAll(".linktext")[0];
		var tpathIndex=tpathList.map(function(d){return paths.indexOf(d[0]);});
		var ttextpath=tpathIndex.map(function(d){return textpaths[d];});
		for (i in ttextpath){
			//var dei=deList[i];
			var dei=tpathList[i][1];
			var tpath=ttextpath[i];
			d3.select(tpath).text("DE gene top "+((dei*100).toFixed(1))+"%")
			.attr("fill","#fff");
		}
		* */
		if (GL.indexOf(deinput)!=-1){		
			plottf(deinput);
		}
		//alert("explore de");
	}

	function explorede_edges(){
		resetPath();
		//alert("this");
		var deinput=d3.select("#deNameEdge").property("value").toUpperCase();

		if (GL.indexOf(deinput)!=-1){
			plottf_edges(deinput);
		}
		//alert("explore de");
	}

	//plot average expression of a TF's targets
	var exploretftarget=function(){
		resetPath(); // ?
		var tfinput=d3.select("#tfName").property("value").toUpperCase();
		var tftargets=regTargetMap[tfinput];

		//get the list of target genes index that are also expressed
		var targetexList=[];
		for (var target of tftargets){
			if (GL.indexOf(target)!=-1) {
				targetexList.push(GL.indexOf(target));
			}
		}

		//for each edge, get the average expression level of expressed target genes
		var avt=[];
		var edgeIDs=[];

		for (var edge of edges){
			var avi=[];
			for (var targetex of targetexList){
				avi.push(edge.expression[targetex]);
			}

			//calculate average over all targets
			var sumi=avi.reduce(function(a,b){return a+b;},0);
			sumi=(sumi/avi.length);

			avt.push(sumi.toFixed(3));
			edgeIDs.push('P'+edge.id);
		}
		barplot(edgeIDs, avt, "Avg. Expression of "+tfinput+" targets");
	}

	// helper function: build target map
	var buildTargetMap=function(geneRegMap){

		var regTargetMap={};
		for (var tf of geneRegMap) {
			var tfname=tf.id;
			regTargetMap[tfname]=tf.targets;
		}

		// for (var i in geneRegMap){
		// 	var gene=exgene[i];
		// 	for (var gt of geneRegMap[i]){
		// 		var gtname=reg[gt].toUpperCase();
		// 		if (!(gtname in regTargetMap)){
		// 			regTargetMap[gtname]=[gene];
		// 		}else{
		// 			regTargetMap[gtname].push(gene);
		// 		}
		// 	}
		// }

		return regTargetMap;
	}

	//reset path
	function resetPath(){
		d3.select("svg").selectAll("path")
		.attr("stroke","#fff");
		
		d3.selectAll(".linktext")
		.text("");
	}

	//reset color bar
	function resetcolorbar(){
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
		
		d3.select(".colorbar").selectAll("text")
		.data([-4,4])
		.text(function(d){return d;});
	}
	
	function setgeneslider(newValue){
		document.getElementById("cellGeneSlider").innerHTML=newValue;
		ptest();
	}

	//zoom slider bar action
	function zoom(newValue){
		document.getElementById("zoomslider").innerHTML=newValue;
		var wd=1000;
		var ht=1200;
		var sv=50;
		var zx=newValue/sv;
		var newwd=wd*zx;
		var newht=ht*zx;
		
		d3.select("#div_svg").select("svg")
		.attr("width",newwd)
		.attr("height",newht)
		.attr("preserveAspectRatio", "none");
	}
	
	//
	function resizecell(){
		var cellsize=document.getElementById("cellsize").value;
		d3.select("svg")
		.selectAll(".cells")
		.attr("r",cellsize);
		
	}

	// show FateBias from FateID
	var showfate=function(){
		var cmat=Object.values(entropy);
		cmat=cmat.filter(function(d){
				if (d>0){
					return true;
				}return false;
			})
		var cmax=Math.max.apply(Math,cmat);
		var cmin=Math.min.apply(Math,cmat);
		
		d3.select("svg").selectAll("g")
		.selectAll(".cells")
		.style("fill",function(d,i){
			var cellcolor=entropy[d.id];
			var cellcolor=(cellcolor-cmin)/(cmax-cmin);
			return hsl_col_perc(cellcolor,0,60);
			})
		.style("opacity",0.3);
	}
	

	//plot the cells
	var drawCells=function(cells){
		
		var errors=cells.map(function(d){return Math.abs(d.location[2]);});
		errors.sort()
		merror=errors[Math.trunc(errors.length/2)];

		var cellabels=[];
		for (var ci of cells){
			if (cellabels.indexOf(ci.type)==-1){
				cellabels.push(ci.type);
			}

		}
		
		var usedcolors=[];
		for (var li in cellabels){
			var lii=cellabels[li];
			var clii=Colors.names[Object.keys(Colors.names)[li]];
			usedcolors.push([lii,clii]);
		}
		
		
		
		circles=d3.select("svg").selectAll("g")
		.selectAll("circle")[0];
		circles=circles.map(function(d){
				return d.__data__;
			})
		
		
		
		var cellsize=3;
		var RSize=16;
		
		d3.select("svg")
		.append("g")
		.selectAll("circle")
		.data(cells)
		.enter()
		.append("circle")
		.attr("class","cells")
		.attr("cx",function(d){
			var ei=edges[d.location[0]];
			var y1=circles[ei.source_node].y;
			var y2=circles[ei.end_node].y;
			var x1=circles[ei.source_node].x;
			var x2=circles[ei.end_node].x;
			var tt=d.location[1];
			var xd=y1+(y2-y1) *tt
			var yd=x1+(x2-x1) *tt;
			var theta=Math.atan((yd-y1)/(xd-x1));
			var R=d.location[2];
			var xr=xd+R*Math.sin(theta);
			var yr=yd+R*Math.cos(theta);
			return xr;
			
			})
		.attr("cy",function(d){
			var ei=edges[d.location[0]];
			var y1=circles[ei.source_node].y;
			var y2=circles[ei.end_node].y;
			var x1=circles[ei.source_node].x;
			var x2=circles[ei.end_node].x;
			var tt=d.location[1];
			var xd=y1+(y2-y1) *tt
			var yd=x1+(x2-x1) *tt;
			var theta=Math.atan((yd-y1)/(xd-x1));
			var cellError=d.location[2];
			var cellError=cellError/merror;
			var R=cellError*RSize;
			var xr=xd+R*Math.sin(theta);
			var yr=yd+R*Math.cos(theta);
			return yr;
			
			})
		.attr("r",cellsize)
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
		.style("fill",function(d){
			//console.log(d.type);
			var ci=cellabels.indexOf(d.type);
			var ci=ci%Object.keys(Colors.names).length;
			var ckey=Object.keys(Colors.names)[ci];
			var cl=Colors.names[ckey];
			/**
			if (d.type=="AT1"){
				return "red";
			}else{
				return "white";
			}
			**/
			
			return cl;
			})
		.style("opacity",0.4);
		
		
		// add the color legend
		d3.select("svg")
		.append("g")
		.attr("class","colorlegend")
		.selectAll("circle")
		.data(usedcolors)
		.enter()
		.append("rect")
		.attr("x",function(d,i){
			return 20;
		})
		.attr("y",function(d,i){
			return 60+i*20;
		})
		.attr("width",20)
		.attr("height",20)
		.style("fill",function(d){
			return d[1];
		})
		
		//add the color text
		d3.select("svg")
		.select(".colorlegend")
		.append("g")
		.selectAll("text")
		.data(usedcolors)
		.enter()
		.append("text")
		.attr("x",function(d,i){
				return 20+20;
			})
		.attr("y",function(d,i){
			return 60+i*20+15;
		})
		.text(function(d){
				return d[0];
		})
		.attr("fill","white");	
		
				
	}
	//draw edge text
	function drawEdgeTexts(TFs){
		//creat TF display svg text 
		createEdgeTF();
		createEdgeETF();
		hideTF();
		hideETF();
		//create DE display SVG text
		createDE();
		//create colorbar 
		createColorBar();
	};
	
	// plot the SVG figure
	function drawTree(root){
		margin = {top: 100, right: 100, bottom: 100, left: 100},  // 600 for plotting for manuscript
		width = 1400 - margin.right - margin.left,
		height =1200 - margin.top - margin.bottom;
		var i = 0;
		var tree = d3.layout.tree();
		tree.size([height,width]);
		tree.separation(function(a,b){
			return (a.parent==b.parent ?1:0.4)/a.depth;
			});
		
		//default bgcolor
		var default_bgcolor="#fff";
		var nodecolor="#a9a9a9"
		var textcolor="#000"
		
		
		var bgcolor;
		var input_bgcolor=d3.select("#bgcolor").value;
		var linecolor="#aaa"
		
		colorisOK  = /(^#[0-9A-F]{6}$)|(^#[0-9A-F]{3}$)/i.test(input_bgcolor);
		
		if (colorisOK){
			bgcolor=input_bgcolor;
		}else{
			bgcolor=default_bgcolor;
		}
		svg = d3.select("#div_svg").append("svg")
			.attr("width", width + margin.right + margin.left)
			.attr("height", height + margin.top + margin.bottom)
			.attr("id","svg")
			.style("background",bgcolor)
			.attr("viewBox","0 0 1400 1200")
			.attr("preserveAspectRatio", "none")
			.append("g")
			.attr("transform", "translate(" + margin.left + "," + margin.top + ")");
			
		var allnodes=tree.nodes(root);
		var links=tree.links(allnodes);
		 
		// Declare the nodes
		var node = svg.selectAll("g.node")
		  .data(nodes)
		  .data(nodes, function(d) { return d.id; });

		// Enter the nodes.
		var nodeEnter = node.enter().append("g")
		  .attr("class", "node")
		  .style("font","8px Arial")
		  .attr("transform", function(d) { 	 
			  return "translate(" + d.y + "," + d.x + ")"; })

		nodeEnter.append("circle")
		  .attr("r", 20)
		  .attr("fill",nodecolor)
		  .attr("stroke","steelblue")
		  .attr("stroke-width","3px")
		  .on("mouseover",inNode)  
		  .on("mouseout",outNode)
		  .on("click",testClick);

		var textEnter=nodeEnter.append("text")
		  .attr("dy", "-2em")
		  .attr("x",-10)
		  .attr("class","nodetext")
		  .attr("fill",textcolor)
		  .style("font","12px Arial")
		  .attr("text-anchor", function(d) { 
			  return  "start"; })
		  .style("fill-opacity", 1)
		  .text(function(d){
			  return 'N'+d.id;});

		//textEnter.call(addwraptext);
		  
		
		
		var link=svg.selectAll("path.link")
				.data(links)
				.enter()
				.append("line")
				.attr("id",function(d,i){return "s"+i;})
				.attr("class","link")
				.attr("stroke","#aaa")
				.attr("stroke-width","2px")
				.attr("x1",function(d){ return d.source.y;})
				.attr("y1",function(d){ return d.source.x;})
				.attr("x2",function(d){ return d.target.y;})
				.attr("y2",function(d){ return d.target.x;});
		
				
		// text along the path
		for (var i in links){
			var ipath=links[i];
			var pid=ipath.target.id-1;
			svg.append("text")
			.attr("x",(ipath.source.y+ipath.target.y)/2)
			.attr("y",(ipath.source.x+ipath.target.x)/2-30)
			.text('P'+pid);
			// 	.style("font-size","20px")
			// .style("font-weight","bold");
		}
		
		
		
		// Define the div for the tooltip
		tooltipsvg = d3.select("body").append("div")	
			.attr("class", "tooltip")
			.append("svg");
			
	    
	    
	}
	//----------------------------------------------------------------------
	//show/hide controls

	//show/hide TF/GE

	function showhideTF(checked,id){
		
		if (checked){
			//pre-check
			var genecheck=d3.select("#div_config").select("#"+id).property("checked");
			if (id=="tfcheck"){
				showTF();
			}
			if (id=="etfcheck"){
				showETF();
			}
		}else{
			if (id=="tfcheck"){
				hideTF();
			}
			if (id=="etfcheck"){
				hideETF();
			}
		}
	}
	
	//----------------------------------------------------------------------
	//hide controls

	//hide TF
	function hideTF(){
		//hide TF text
		d3.select("svg").selectAll(".TFText")
		.style("opacity",0);
		
		d3.select("svg").selectAll(".TFText")
		.selectAll("tspan")
		.style("font-size",0)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		//hide TF color bar
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
	}
	
	function hideETF(){
		//hide TF text
		d3.select("svg").selectAll(".ETFText")
		.style("opacity",0);
		
		d3.select("svg").selectAll(".ETFText")
		.selectAll("tspan")
		.style("font-size",0)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		//hide TF color bar
		d3.select("svg").select(".colorbar")
		.style("opacity",0);
	}
	//---------------------------------------------------------------------
	//show controls

	//show TF 
	function showTF(){
		
		//change the opacity
		d3.select("svg").selectAll(".TFText")
		.style("opacity",1);
		//change the text size
		d3.select("svg").selectAll(".TFText")
		.selectAll("tspan")
		.style("font-size",8)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		
		resetcolorbar();
		d3.select("svg").select(".colorbar")
		.style("opacity",1);
		
		
	}
	
	function showETF(){
		
		//change the opacity
		d3.select("svg").selectAll(".ETFText")
		.style("opacity",1);
		//change the text size
		d3.select("svg").selectAll(".ETFText")
		.selectAll("tspan")
		.style("font-size",8)
		.attr("x","-8em")
		.attr("dy","-1em");
		
		
		resetcolorbar();
		d3.select("svg").select(".colorbar")
		.style("opacity",1);
		
	}

	//----------------------------------------------------------------------
	//create svg text for showing DE gens

	function createDE(){
		var gc_enter=svg.selectAll("g")
		.append("text")
		.attr("class","DEText")
		.attr("opacity",0)
		.attr("x","-10em")
		.attr("y","-3em")
		gc_enter.call(addDEText,20);
	}

	//add DE texts 
	function addDEText(text,showcutoff){
		text.each(function(d){
			var chosenEdge=null;
			for (var edge of edges){
				if (nodes[edge.to]==d){
					chosenEdge=JSON.parse(JSON.stringify(edge));
				}
			}
			if (chosenEdge!=null){
				var de=chosenEdge.de;
				de=de.map(function(x){return x.toUpperCase();});
				var de=de.splice(0,20).reverse();
				for (var gd of de){
					var pnodeE=nodes[chosenEdge.from].E;
					var cnodeE=d.E;
					var tfindex=GL.indexOf(gd);
					var fc=cnodeE[tfindex]-pnodeE[tfindex];
					fc=fc.toFixed(2);
					nfc=Math.min(Math.max(fc,-4),4);
					var per=(4-nfc)/(4+4);
					var cl=hsl_col_perc(per,0,120);
					d3.select(this).append("tspan")
					.style("font-size",8)
					.attr("x","-8em")
					.attr("dy","-1em")
					.style("fill",cl)
					.text(gd+' ('+fc+')');
				}
			}
		});
	}
	
	
	var drag=d3.behavior.drag()
			.on("drag",dragmove);
			
	var drage=d3.behavior.drag()
			.on("drag",edragmove);
			
	function dragmove(d){
		
		var dragx = d3.event.x;
		var dragy = d3.event.y;
		
		var tx=d3.select(this).select("text").node().x.baseVal[0].value;
		var ty=d3.select(this).select("text").node().y.baseVal[0].value;
		
		var offsetx=dragx-tx;
		var offsety=dragy-ty;
		
		d3.select(this).attr("transform","translate(" + offsetx+","+ offsety +")");
		
		d=d[0];
		var ei=edges[d.location[0]];
		var y1=circles[ei.source_node].y;
		var y2=circles[ei.end_node].y;
		var x1=circles[ei.source_node].x;
		var x2=circles[ei.end_node].x;
		var tt=d.location[1];
		var xd=y1+(y2-y1) *tt+margin.left;
		var yd=x1+(x2-x1) *tt+margin.top;
		
		d3.select("svg")
		.select("line#line"+this.id)
		.style("stroke","#aaa")
		.attr("class","TFText")
		.attr("x1",xd)
		.attr("y1",yd)
		.attr("x2",dragx+margin.left)
		.attr("y2",dragy+margin.top);
		
	}
	
	
	function edragmove(d){
		
		var dragx = d3.event.x;
		var dragy = d3.event.y;
		
		var tx=d3.select(this).select("text").node().x.baseVal[0].value;
		var ty=d3.select(this).select("text").node().y.baseVal[0].value;
		
		var offsetx=dragx-tx;
		var offsety=dragy-ty;
		
		d3.select(this).attr("transform","translate(" + offsetx+","+ offsety +")");
		
		d=d[0];
		var ei=edges[d.location];
		var y1=circles[ei.source_node].y;
		var y2=circles[ei.end_node].y;
		var x1=circles[ei.source_node].x;
		var x2=circles[ei.end_node].x;
		var tt=1;
		var xd=y1+(y2-y1) *tt+margin.left;
		var yd=x1+(x2-x1) *tt+margin.top;
		
		d3.select("svg")
		.select("line#eline"+this.id.slice(1))
		.attr("class","ETFText")
		.style("stroke","#aaa")
		.attr("x1",xd)
		.attr("y1",yd)
		.attr("x2",dragx+margin.left)
		.attr("y2",dragy+margin.top);
		
	}

	
	function createEdgeTF(){
		var usedCors={};
		var TFGroups={};
		for (var tf of TFs){
			ei=edges[tf.location[0]];
			tt=tf.location;
			var key=ei.source_node+","+ei.end_node+","+tf.location[1];
			if (!(key in TFGroups)){
				TFGroups[key]=[tf];
			}else{
				TFGroups[key].push(tf);
			}
		}
		var TFGroupList=Object.values(TFGroups);
		
		var gid=0;
		d3.select("svg").append("g")
		.selectAll("g")
		.data(TFGroupList)
		.enter()
		.append("g")
		.call(drag)
		.attr("class","tfs")
		.attr("id",function(d,i){
				return d[0].id;
			});
		
		d3.select("svg")
		.append("g")
		.selectAll("line")
		.data(TFGroupList)
		.enter()
		.append("line")
		.attr("id",function(d,i){
			//console.log(d[0].id);
			return "line"+d[0].id;
			
			})
		
		
		d3.selectAll("g.tfs")
		.each(function(d,i){
			d3.select(this)
			.selectAll("text")
			.data(d)
			.enter()
			.append("text")
			.attr("class","TFText")
			.attr("x",function(d,i){
				var ei=edges[d.location[0]];
				var y1=circles[ei.source_node].y;
				var y2=circles[ei.end_node].y;
				var x1=circles[ei.source_node].x;
				var x2=circles[ei.end_node].x;
				var tt=d.location[1];
				var xd=y1+(y2-y1) *tt
				var yd=x1+(x2-x1) *tt;
				return xd;
			})
			.attr("y",function(d,i){
				var ei=edges[d.location[0]];
				var y1=circles[ei.source_node].y;
				var y2=circles[ei.end_node].y;
				var x1=circles[ei.source_node].x;
				var x2=circles[ei.end_node].x;
				var tt=d.location[1];
				var xd=y1+(y2-y1) *tt
				var yd=x1+(x2-x1) *tt;
				var offset=25
				if (!([xd,yd] in usedCors)){
					usedCors[[xd,yd]]=0;
					return yd+offset;
				}else{
					usedCors[[xd,yd]]+=1;
					yd+=16*usedCors[[xd,yd]];
					return yd+offset;
				}
				
				
			})
			.text(function(d,i){
				 return 'path ' + d.location[0] +" time "+ d.location[1] + ':' + d.tf_name.toUpperCase()+"("+d["p-value"].toExponential(2)+")";
				 //console.log(d);
			})
			.style("fill","white")
			.style("font-size","10px")
			.style("font-weight","bold")
			.style("text-anchor","start")
			.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
			});

		
	}	
			
	//create svg text for showing TFs
	function createEdgeTF1(){
		var usedCors={};
		d3.select("svg").append("g")
		.attr("opacity",0.5)
		.attr("class","tfs")
		.selectAll("text")
		.data(TFs)
		.enter()
		.append("text")
		.attr("x",function(d,i){
			var ei=edges[d.location[0]];
			var y1=circles[ei.source_node].y;
			var y2=circles[ei.end_node].y;
			var x1=circles[ei.source_node].x;
			var x2=circles[ei.end_node].x;
			var tt=d.location[1];
			var xd=y1+(y2-y1) *tt
			var yd=x1+(x2-x1) *tt;
			return xd;
		})
		.attr("y",function(d,i){
			var ei=edges[d.location[0]];
			var y1=circles[ei.source_node].y;
			var y2=circles[ei.end_node].y;
			var x1=circles[ei.source_node].x;
			var x2=circles[ei.end_node].x;
			var tt=d.location[1];
			var xd=y1+(y2-y1) *tt
			var yd=x1+(x2-x1) *tt;
			var offset=25
			if (!([xd,yd] in usedCors)){
				usedCors[[xd,yd]]=0;
				return yd+offset;
			}else{
				usedCors[[xd,yd]]+=1;
				yd+=16*usedCors[[xd,yd]];
				return yd+offset;
			}
			
			
		})
		.text(function(d,i){
			 return d.tf_name;
			 //console.log(d);
		})
		.style("fill","white")
		.style("font-size","16px")
		.style("font-weight","bold")
		.style("text-anchor","start")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
		.call(drag);
		
		
	}
	
	function createEdgeETF(){
		
		var textcolor="#000000"
		var usedCors={};
		var TFGroups={};
		for (var tf of ETFs){
			ei=edges[tf.location];
			tt=tf.location;
			var key=ei.source_node+","+ei.end_node+","+'1';
			if (!(key in TFGroups)){
				TFGroups[key]=[tf];
			}else{
				TFGroups[key].push(tf);
			}
		}
		var TFGroupList=Object.values(TFGroups);
		
		var gid=0;
		d3.select("svg").append("g")
		.selectAll("g")
		.data(TFGroupList)
		.enter()
		.append("g")
		.call(drage)
		.attr("class","tfs")
		.attr("id",function(d,i){
				return 'e'+d[0].id;
			});
		
		d3.select("svg")
		.append("g")
		.selectAll("line")
		.data(TFGroupList)
		.enter()
		.append("line")
		.attr("id",function(d,i){
			//console.log(d[0].id);
			return "eline"+d[0].id;
			
			})
		
		
		d3.selectAll("g.tfs")
		.each(function(d,i){
			d3.select(this)
			.selectAll("text")
			.data(d)
			.enter()
			.append("text")
			//.attr("class","etfs")
			.attr("class","ETFText")
			.attr("x",function(d,i){
				var ei=edges[d.location];
				var y1=circles[ei.source_node].y;
				var y2=circles[ei.end_node].y;
				var x1=circles[ei.source_node].x;
				var x2=circles[ei.end_node].x;
				var tt=1;
				var xd=y1+(y2-y1) *tt
				var yd=x1+(x2-x1) *tt;
				return xd;
			})
			.attr("y",function(d,i){
				var ei=edges[d.location];
				var y1=circles[ei.source_node].y;
				var y2=circles[ei.end_node].y;
				var x1=circles[ei.source_node].x;
				var x2=circles[ei.end_node].x;
				var tt=1;
				var xd=y1+(y2-y1) *tt
				var yd=x1+(x2-x1) *tt;
				var offset=25
				if (!([xd,yd] in usedCors)){
					usedCors[[xd,yd]]=0;
					return yd+offset;
				}else{
					usedCors[[xd,yd]]+=1;
					yd+=16*usedCors[[xd,yd]];
					return yd+offset;
				}
				
				
			})
			.text(function(d,i){

				 return 'eTF: path ' + d.location + ':' +d.etf_name.toUpperCase()+"("+d["p-value"].toExponential(2)+")";
				 //console.log(d);
			})
			.style("fill",textcolor)
			.style("font-size","10px")
			.style("font-weight","bold")
			.style("text-anchor","start")
			.attr("transform", "translate(" + margin.left + "," + margin.top + ")")
			});

		
	}	
	
	// add addTFText
	function addTFText(){
		
	}


	// create color bar 
	function createColorBar(){
		var textcolor="#000000"
		var cx=[];
		for (var ci=1;ci>0;ci-=0.005){
			cx.push(ci);
		}
		var xstart=10;
		var ystart=15;
		var xspan=1;
		var xspan_height=10;
		var xspan_width=4;
		d3.select("svg").append("g")
		.attr("opacity",0)
		.attr("class","colorbar")
		.selectAll("rect")
		.data(cx)
		.enter()
		.append("rect")
		.attr("x",function(d,i){
			var xx=xstart+i*xspan;
			return xx;
			})
		.attr("y",ystart)
		.attr("width",xspan_width)
		.attr("height",xspan_height)
		.attr("fill",function(d){
			var cl=hsl_col_perc(d,0,120);
			return cl;
			});
		xmm=[0,1]
		d3.select("svg").select(".colorbar").append("g")
		.selectAll("text")
		.data(xmm)
		.enter()
		.append("text")
		.style("fill",textcolor)
		.attr("x",function(d,i){
				xx=xstart+d*xspan*cx.length;
				return xx;
			})
		.attr("y",ystart-2)
		.text(function(d){
			if(d==0){
				return -4;
			}else{
				return 4;
				}})
		.attr("text-anchor",function(d,i){
				if(d==0){
					return "start";
				}else{
					return "end";
				}
			});
	}

	//color convertion
	function hsl_col_perc(percent,start,end) {

		 var a=percent;
		 b = end*a;
		 c = b+start;

		//Return a CSS HSL string
		return 'hsl('+c+',100%,50%)';
	}

	//get the expression for given node
	function getNodeEx(node){
		var ne=new Array(cells[node.CELL[0]].E.length);
		ne.fill(0);
		
		for(var cell of node.CELL){
			var ce=cells[cell].E;
			for (var x in ce){
				ne[x]+=ce[x];
			}
		}
		for (x in ne){
			ne[x]/=node.CELL.length;
		}
		return ne;
	}

	// wrap multi- line text 
	function addwraptext(text){
		text.each(function(d){
			CELL=d.CELL;
			var ctdict={};
			for(x of CELL){
				xlabel=cells[x].typeLabel;
				if (xlabel in ctdict==false){
					ctdict[xlabel]=1;
				}else{
					ctdict[xlabel]+=1;
				}
			}
			d3.select(this).append("tspan")
				.attr("x",20)
				.attr("dy","1em")
				.text(function(d){
					return "ID: E"+d.T+"_"+d.ID;
					})
			// add diff stage
			d3.select(this).append("tspan")
				.attr("x",20)
				.attr("dy","1em")
				.text(function(d){
					return "Differentiation: "+(d.D*100).toFixed(1)+"%";
				})
				
			for ( var key in ctdict){
				d3.select(this).append("tspan")
				.attr("x",function(d){ return 20;})
				.attr("dy","1em")
				.text(key+","+ctdict[key]);
			}
		});
	}


	// get Node text
	function getNodeText(d){
		return d.ID;
	}

	//node click function
	function testClick(){
		var newW = open('','_blank','height=600,width=800,left=200,top=200,scrollbars=yes')
		newW.document.write("<head><title>Edge TF/Gene</title></head><body></body>");
		newW.document.title="loading...";

		//adding TF details
		var tdiv=d3.select(newW.document.body)
		.style("background","white")
		.append("div")
		.style("padding-left","50px")
		.style("padding-top","50px")
		.attr("width",400)
		.attr("height",600)
		.attr("class","div_table_edge");
		
		
		cnode=this.__data__;
		// cell details at selected node 
		var CellIDs=[["ID","Label"]];
		for (var icell of cnode.CELL){
			var cid=cells[icell].ID;
			var clabel=cells[icell].typeLabel;
			CellIDs.push([cid,clabel]);
		}
		
		tdiv
		.append("p")
		.text("Table 0. Cells at selected node: "+cnode.T+"_"+cnode.ID);
		
		
		createTable(tdiv,"celltable",CellIDs);
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"celltable","celldlink");
			});
			
		tdiv.append("a")
		.attr("href","#")
		.attr("id","celldlink");
		
		
		//TF details for edge ending at selected node
		pnode=cnode.parent;
		var chosenEdge;
		for (var edge of edges){
			if (nodes[edge.to]==cnode){
				chosenEdge=JSON.parse(JSON.stringify(edge));
			}
		}
		
		var resList=[]; //TF details
		resList.push(["TF","p-value","TF fold change","Mean target fold change","Mean DE target fold change","edgeID"]);

		resList=exportTFEdge(resList,chosenEdge);
		
		tdiv
		.append("p")
		.text("Table 1. TF details for the edge ending at selected node:"+cnode.T+"_"+cnode.ID);
		
		
		createTable(tdiv,"tftable",resList);
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"tftable","tfdlink");
			});
			
		tdiv.append("a")
		.attr("href","#")
		.attr("id","tfdlink");
		
		
		////TF details for path ending at selected node
		cnode=this.__data__;
		pnode=cnode.parent;
		var chosenPath=[]
		
		while(pnode!="null"){
			var ec=cnode.endEdge;
			cnode=pnode;
			pnode=pnode.parent;
			chosenPath.push(ec);
		}
		

		var AllTFList=[];
		AllTFList.push(["TF","p-value","TF fold change","Mean target fold change","Mean DE target fold change","edgeID"]);

		for (var cedge of chosenPath){
			AllTFList=exportTFEdge(AllTFList,JSON.parse(JSON.stringify(cedge)));
			
		}

		AllTFList.sort(function(a,b){return a[1]-b[1];});
		tdiv
		.append("p")
		.text("Table 2. TFs along the path ending at selected node:"+cnode.T+"_"+cnode.ID);
		
		createTable(tdiv,"alltftable",AllTFList);
		
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"alltftable","alltfdlink");
			});
			
		tdiv.append("a")
		.attr("href","#")
		.attr("id","alltfdlink")
		.text("");
		
		//adding DE details
		cnode=this.__data__;
		pnode=cnode.parent;
		var resDEList=[]; //DE details
		resDEList.push(["DE gene","Expression_E"+pnode.T+"_"+pnode.ID,"Expression_"+cnode.T+"_"+cnode.ID,"Fold change","edgeID"]);

		var resDEList=exportDEEdge(resDEList,chosenEdge);
		tdiv.append("p")
		.text("Table 3. DE gene details for edge ending at selected node:"+cnode.T+"_"+cnode.ID)
		
		createTable(tdiv,"detable",resDEList);
		tdiv.append("button")
		.text("click to create the excel table")
		.on("click",function(){
				table2XLS(newW,"detable","dedlink");
			});
		tdiv.append("a")
		.attr("href","#")
		.attr("id","dedlink");
		
		//loading complete
		newW.document.title="loading complete";
		
	}

	//function export tf details of a edge
	function exportTFEdge(resList,chosenEdge){
		var tf,tfID,pv,rank,fc,tfc;
		var etf=chosenEdge.etf;
		var de=chosenEdge.de;
		de=de.map(function(d){return d.toUpperCase();});
		cnode=nodes[chosenEdge.to];
		pnode=nodes[chosenEdge.from];
		var edgeid='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;
		for (i in etf){
			tf=etf[i];
			var tfID=tf[1].toUpperCase();
			var tfp=tf[0].toExponential(2);
			var tfrank=i;
			var fc=cnode.fc[tfID].toFixed(3);
			var target=dTD[tfID].filter(function(d){return d in cnode.fc;});
			var tfc=target.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
			tfc=(tfc.reduce(function(a,b){return a+b;})/tfc.length).toFixed(3);
			var dtfc=target.filter(function(d){
					if (de.indexOf(d)!=-1){
						return true;
					}else{
						return false;
					}
				});
			dtfc=dtfc.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
			dtfc=(dtfc.reduce(function(a,b){return a+b;})/dtfc.length).toFixed(3);
			var tfresList=[tfID,tfp,fc,tfc,dtfc,edgeid];
			resList.push(tfresList);
		}
		return resList;
	}

	function exportDEEdge(resDEList,chosenEdge){
		var cnode=nodes[chosenEdge.to];
		var pnode=nodes[chosenEdge.from];
		
		var edgeID='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;


		var etf=chosenEdge.etf;
		var de=chosenEdge.de.map(function(d){return d.toUpperCase();});
		
		for (i in de){
			g=de[i];
			ef=pnode.E[GL.indexOf(g)];
			et=cnode.E[GL.indexOf(g)];
			fc=et-ef;
			resDEList.push([g,ef.toFixed(3),et.toFixed(3),fc.toFixed(3),edgeID]);
		}
		return resDEList;
	}

	//append a table to the cant
	function createTable(cant,tableid,data){
		cant.append("table")
		.attr("id",tableid)
		.style("border","1px solid")
		.style("border-collapse","collapse")
		.selectAll("tr")
		.data(data)
		.enter()
		.append("tr")
		.style("border","1px solid")
		.selectAll("td")
		.data(function(d){return d;})
		.enter()
		.append("td")
		.style("border","1px solid")
		.text(function(d){return d;});
	}


	//action on mouse out nodeID: E1_16.0_0Proliferative AT2 Early Precursor,6Proliferative Bi-potential Precursor,1

	function outNode(){
		tooltipsvg
		.style("opacity", 0)
		.attr("width",0)
		.attr("height",0)
		.selectAll("*").remove();;
	}

	function inNode(){
		if (d3.select("#tooltipcheck").property("checked")==false){
			return false;
		}
		var CELL=this.__data__.CELL;
		var ctdict={}
		for(x of CELL){
			xlabel=cells[x].typeLabel;
			if (xlabel in ctdict==false){
				ctdict[xlabel]=1;
			}else{
				ctdict[xlabel]+=1;
			}
		}

		var ltdict=[]
		var asum=0;
		for (key in ctdict){
			asum+=ctdict[key];
			ltdict.push([key,ctdict[key]])
		}
		
		var cellCounter=0;
		
		var arc=d3.svg.arc()
				.innerRadius(100)
				.outerRadius(150)
				.startAngle(function(d,i){
					var sa=2*Math.PI*(cellCounter/asum);
					return sa;
					})
				.endAngle(function(d,i){
					cellCounter+=d[1];
					var ea=2*Math.PI*(cellCounter/asum);
					return ea;
					});
			
		tooltipsvg
			.style("left",(d3.event.pageX+30)+"px")
			.style("top",(d3.event.pageY-28)+"px")
			.style("position","absolute")
			.style("opacity",1)
			.style("background","gray")
			//.style("border","1px solid")
			.attr("class","tooltipsvg")
			.attr("width",300)
			.attr("height",300)
		
		var colorCodes=["red","blue","purple","green","lime","magenta","orange","olive","cyan","yellow"];
		var colorDex=0;
		tooltipsvg.selectAll("path")
			.data(ltdict)
			.enter()
			//.append("circle")
			//.attr("r",10);
			.append("path")
			.attr("d",arc)
			.attr("transform","translate (150,150)")
			.attr("fill",function(d,i) { 
				colorDex=colorDex % 10;
				var colorChoosen=colorCodes[colorDex];
				colorDex+=1;
				return colorChoosen;}
				)
				
		cellCounter=0;
		colorDex=0;
		tooltipsvg.selectAll("text")
			.data(ltdict)
			.enter()
			.append("text")
			.style("border","1px solid")
			.attr("transform",function(d,i){
				var tx=60;
				var ty=100+i*20;
				return "translate ("+tx+","+ty+")";
				})
			.attr("dy",20)
			.text(function(d,i){return d[0]+':'+d[1];})
			.style("font-size","10px")
			.style("fill",function(d,i) { 
				colorDex=colorDex % 10;
				var colorChoosen=colorCodes[colorDex];
				colorDex+=1;
				return colorChoosen;}
				);
	}



	//action on mouse over node
	function oldnodemouseoverfunction(){
		
		var CELL=this.__data__.CELL;
		ctdict={}
		for(x of CELL){
			xlabel=cells[x].typeLabel;
			if (xlabel in ctdict==false){
				ctdict[xlabel]=1;
			}else{
				ctdict[xlabel]+=1;
			}
		}
		
		tooltipdiv
		.html(JSON.stringify(ctdict))
		.style("left",(d3.event.pageX)+"px")
		.style("top",(d3.event.pageY-28)+"px")
		.style("position","absolute")
		.style("opacity",0.9)
		.transition()
		.duration(200);
		
	}

	function nodeClick(){
		var newW=window.open("","","width=300,height=250");
		newW.document.write("<html><body></body></html>");
		var svg = d3.select(newW.document.body)
			.append("svg")
			.attr("width",400)
			.attr("height",200)
			.append("circle")
			.attr("cx",200)
			.attr("cy",100)
			.attr("r",20)
	}

	//create tree structure
	function buildTree(nodes,edges){
		for (var nodeIndex in nodes){
			var childrenList=[];
			for (var edgeIndex in edges){
				if (edges[edgeIndex].source_node==nodeIndex){
					childrenList.push(nodes[edges[edgeIndex].end_node]);
				}
				
			}
			nodes[nodeIndex].name=nodes[nodeIndex].id;
			nodes[nodeIndex].children=childrenList;
		}
		
		return nodes;
		
	}

	// parse json file
	function parseJSON(data){
		// note here the names are not the same as in onload ; for example, dTD may correspond to TFs
		data=JSON.parse(data);
		var GL=data[0];
		var CellList=data[1];
		var NodeList=data[2];
		var EdgeList=data[3];
		var dTD=data[4];	
		return [GL,CellList,NodeList,EdgeList,dTD,data[5]];
	};


	//---------------------------------------------------------------------
	//download functions here

	//download tf
	function downloadde(){
		plswait("dedownloadlink");
		window.setTimeout(creatededownload,10);
	}

	function creatededownload(){
		var resList=[];
		resList.push(["DE gene","Expression_from","Expression_to","Fold change","edgeID"]);

		for (var edge of edges){
			var etf=edge.etf;
			var de=edge.de;
			de=de.map(function(d){return d.toUpperCase();});
			var cnode=nodes[edge.to];
			var pnode=nodes[edge.from];
			var edgeid='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;
			for (var g of de){
				var gindex=GL.indexOf(g);
				var efrom=pnode.E[gindex];
				var eto=pnode.E[gindex];
				var fc=eto-efrom;
				var deresList=[g,efrom,eto,fc,edgeid];
				resList.push(deresList);
			}
		}
		//resList=resList.splice(0,10000);
		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#dedownloadlink")
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}
	
	function downloadcell(){
		var nw=window.open("","","width=600,height=600");
			var tdiv = d3.select(nw.document.body)
			.append("div");
			var delist1=s1.map(function(d){return [1,d.__data__.id];});
			var delist2=s2.map(function(d){return [2,d.__data__.id];});
			var delist=delist1.concat(delist2);
			createTable(tdiv,"celltable",[['Group_ID','Cell_ID']].concat(delist));
		
	}
	
	//download tf
	function downloadtf(){
		plswait("tfdownloadlink");
		window.setTimeout(createtfdownload2,10);
	}

	//download tf
	function createtfdownload2(){
		var resList=[];
		resList.push(["TF","pathID", "Location on edge","p-value"]);
		for (var tf of TFs){
			var tfID=tf.tf_name.toUpperCase();
			var tfp=tf["p-value"].toExponential(2);
			var tfrank=tf.id;
			var tfEdge=tf.location[0];
			var tfLocation=tf.location[1];
			var tfresList=[tfID,tfEdge,tfLocation, tfp];
			resList.push(tfresList);
		}

		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#tfdownloadlink")
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}

	//download etf
	function downloadetf(){
		plswait("etfdownloadlink");
		window.setTimeout(createetfdownload,10);
	}

	//download etf
	function createetfdownload(){
		var resList=[];
		resList.push(["eTF","pathID", "fold change compared to parent","fold change compared to children","p-value"]);
		for (var tf of ETFs){
			var tfID=tf.etf_name.toUpperCase();
			var tfp=tf["p-value"].toExponential(2);
			var tfFcp=tf.fcp;
			var tfFcs=tf.fcs;
			var tfEdge=tf.location;
			var tfresList=[tfID,tfEdge,tfFcp, tfFcs, tfp];
			resList.push(tfresList);
		}

		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#etfdownloadlink")
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}

	//download tf
	function createtfdownload(){
		var resList=[];
		resList.push(["TF","p-value","TF fold change","Mean target fold change","Mean DE target fold change","edgeID"]);
		for (var edge of edges){
			var etf=edge.etf;
			var de=edge.de;
			de=de.map(function(d){return d.toUpperCase();});
			var cnode=nodes[edge.to];
			var pnode=nodes[edge.from];
			var edgeid='E'+pnode.T+'_'+pnode.ID+'->E'+cnode.T+"_"+cnode.ID;
			for (var i in etf){
				tf=etf[i];
				var tfID=tf[1].toUpperCase();
				var tfp=tf[0].toExponential(2);
				var tfrank=i;
				var fc=cnode.fc[tfID].toFixed(3);
				var target=dTD[tfID].filter(function(d){return d in cnode.fc;});
				var tfc=target.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
				tfc=(tfc.reduce(function(a,b){return a+b;})/tfc.length).toFixed(3);
				var dtfc=target.filter(function(d){
						if (de.indexOf(d)!=-1){
							return true;
						}else{
							return false;
						}
					});
				dtfc=dtfc.map(function(d){return cnode.E[GL.indexOf(d)]-pnode.E[GL.indexOf(d)];});
				dtfc=(dtfc.reduce(function(a,b){return a+b;})/dtfc.length).toFixed(3);
				var tfresList=[tfID,tfp,fc,tfc,dtfc,edgeid];
				resList.push(tfresList);
			}
		}
		var outString=List2TSV(resList);
		var blob=new Blob([outString],{type:'application/vnd.ms-excel'});
		outurl=window.URL.createObjectURL(blob);
		d3.select("#tfdownloadlink")
		.attr("href",outurl)
		.attr("download","download.xls")
		.text("Ready,Click to download");
	}

	function plswait(id){
		document.getElementById(id).innerHTML="wait...";
		d3.select("#"+id)
		.text("Generating file,please wait...");	
	}

	//convert 2D-list to  TSV
	function List2TSV(LST){
		var out="";
		for (var i in LST){
			var LSTI=LST[i].join("	");
			out+=LSTI+'\n';
		}
		return out;
	}
	//write download link


	//download figure

	function downloadfig(){
		svgToCanvas();

	}
	
	// save svg
	function saveSvg(svgEl, name) {
		svgEl.setAttribute("xmlns", "http://www.w3.org/2000/svg");
		var svgData = svgEl.outerHTML;
		var preface = '<?xml version="1.0" standalone="no"?>\r\n';
		var svgBlob = new Blob([preface, svgData], {type:"image/svg+xml;charset=utf-8"});
		var svgUrl = URL.createObjectURL(svgBlob);
		var downloadLink = document.createElement("a");
		downloadLink.href = svgUrl;
		downloadLink.download = name;
		document.body.appendChild(downloadLink);
		downloadLink.click();
		document.body.removeChild(downloadLink);
	}
	
	//svgToCanvas
	function svgToCanvas(){
		//get svg element.
		var svg = document.getElementById("svg");
		saveSvg(svg, 'download.svg')
	}

	//svgToCanvas
	/*
	function svgToCanvas(){
		//get svg element.
		var svg = document.getElementById("svg");
		saveSvg(svg, 'download.svg')
		//get svg source.
		var serializer = new XMLSerializer();
		var source = serializer.serializeToString(svg);

		//add name spaces.
		if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)){
			source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');
		}
		if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)){
			source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');
		}

		//add xml declaration
		source = '<?xml version="1.0" standalone="no"?>\n' + source;

		//convert svg source to URI data scheme.
		var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(source);

		//set url value to a element's href attribute.
		document.getElementById("downloadlink").href = url;
		document.getElementById("downloadlink").innerHTML="Figure ready,Right click me to save!"
		//you can download svg file by right click menu	
	}
	* */

	//convert table to excel download link
	//dlink: link id
	function table2XLS(newW,tlink,dlink){
		var datatype="data:application/vnd.ms-excel";
		var table_div=newW.document.getElementById(tlink).outerHTML;
		var url='data:application/vnd.ms-excel,'+encodeURIComponent(table_div);
		newW.document.getElementById(dlink).href=url;
		newW.document.getElementById(dlink).download="download.xls";
		newW.document.getElementById(dlink).innerHTML="Table ready!click to save! (please add the right file extension .xls if not prompted)"
	}
	
	
	// colors
	function getColors(){
			var colorObject={
          "scolor1": "#4dc3cfff",
          "scolor2": "#ed1d24ff",
          "scolor3": "#f46e16ff",
          "scolor4": "#8f55b1ff",
          "scolor5": "#37a53eff",
          "scolor6": "#3974bbff",
		  "burlywood": "#deb887",
		  "cadetblue": "#5f9ea0",
		  "chartreuse": "#7fff00",
		  "chocolate": "#d2691e",
		  "coral": "#ff7f50",
		  "cornflowerblue": "#6495ed",
		  "cornsilk": "#fff8dc",
		  "crimson": "#dc143c",
		  "darkblue": "#00008b",
		  "darkcyan": "#008b8b",
		  "darkgoldenrod": "#b8860b",
		  "darkgray": "#a9a9a9",
		  "darkgreen": "#006400",
		  "darkgrey": "#a9a9a9",
		  "darkkhaki": "#bdb76b",
		  "darkmagenta": "#8b008b",
		  "darkolivegreen": "#556b2f",
		  "darkorange": "#ff8c00",
		  "darkorchid": "#9932cc",
		  "darkred": "#8b0000",
		  "darksalmon": "#e9967a",
		  "darkseagreen": "#8fbc8f",
		  "darkslateblue": "#483d8b",
		  "darkslategray": "#2f4f4f",
		  "darkslategrey": "#2f4f4f",
		  "darkturquoise": "#00ced1",
		  "darkviolet": "#9400d3",
		  "deeppink": "#ff1493",
		  "deepskyblue": "#00bfff",
		  "dimgray": "#696969",
		  "dimgrey": "#696969",
		  "dodgerblue": "#1e90ff",
		  "red": "#ff0000",
		  "orange": "#ffa500",
          "purple": "#800080",
		  "green": "#008000",
		  "blue": "#0000ff",
		  "yellow": "#ffff00",
          "cyan": "#00ffff",
		  "pink": "#ffc0cb",
		  "black": "#000000",
		  "brown": "#a52a2a",
		  "firebrick": "#b22222",
		  "floralwhite": "#fffaf0",
		  "forestgreen": "#228b22",
		  "fuchsia": "#ff00ff",
		  "gainsboro": "#dcdcdc",
		  "ghostwhite": "#f8f8ff",
		  "gold": "#ffd700",
		  "goldenrod": "#daa520",
		  "gray": "#808080",
		  "greenyellow": "#adff2f",
		  "grey": "#808080",
		  "honeydew": "#f0fff0",
		  "hotpink": "#ff69b4",
		  "indianred": "#cd5c5c",
		  "indigo": "#4b0082",
		  "ivory": "#fffff0",
		  "khaki": "#f0e68c",
		  "lavender": "#e6e6fa",
		  "lavenderblush": "#fff0f5",
		  "lawngreen": "#7cfc00",
		  "lemonchiffon": "#fffacd",
		  "lightblue": "#add8e6",
		  "lightcoral": "#f08080",
		  "lightcyan": "#e0ffff",
		  "lightgoldenrodyellow": "#fafad2",
		  "lightgray": "#d3d3d3",
		  "lightgreen": "#90ee90",
		  "lightgrey": "#d3d3d3",
		  "lightpink": "#ffb6c1",
		  "lightsalmon": "#ffa07a",
		  "lightseagreen": "#20b2aa",
		  "lightskyblue": "#87cefa",
		  "lightslategray": "#778899",
		  "lightslategrey": "#778899",
		  "lightsteelblue": "#b0c4de",
		  "lightyellow": "#ffffe0",
		  "lime": "#00ff00",
		  "limegreen": "#32cd32",
		  "linen": "#faf0e6",
		  "magenta": "#ff00ff",
		  "maroon": "#800000",
		  "mediumaquamarine": "#66cdaa",
		  "mediumblue": "#0000cd",
		  "mediumorchid": "#ba55d3",
		  "mediumpurple": "#9370db",
		  "mediumseagreen": "#3cb371",
		  "mediumslateblue": "#7b68ee",
		  "mediumspringgreen": "#00fa9a",
		  "mediumturquoise": "#48d1cc",
		  "mediumvioletred": "#c71585",
		  "midnightblue": "#191970",
		  "mintcream": "#f5fffa",
		  "mistyrose": "#ffe4e1",
		  "moccasin": "#ffe4b5",
		  "navajowhite": "#ffdead",
		  "navy": "#000080",
		  "oldlace": "#fdf5e6",
		  "olive": "#808000",
		  "olivedrab": "#6b8e23",
		  "orangered": "#ff4500",
		  "orchid": "#da70d6",
		  "palegoldenrod": "#eee8aa",
		  "palegreen": "#98fb98",
		  "paleturquoise": "#afeeee",
		  "palevioletred": "#db7093",
		  "papayawhip": "#ffefd5",
		  "peachpuff": "#ffdab9",
		  "peru": "#cd853f",
		  "plum": "#dda0dd",
		  "powderblue": "#b0e0e6",
		  "rebeccapurple": "#663399",
		  "rosybrown": "#bc8f8f",
		  "royalblue": "#4169e1",
		  "saddlebrown": "#8b4513",
		  "salmon": "#fa8072",
		  "sandybrown": "#f4a460",
		  "seagreen": "#2e8b57",
		  "seashell": "#fff5ee",
		  "sienna": "#a0522d",
		  "silver": "#c0c0c0",
		  "skyblue": "#87ceeb",
		  "slateblue": "#6a5acd",
		  "slategray": "#708090",
		  "slategrey": "#708090",
		  "snow": "#fffafa",
		  "springgreen": "#00ff7f",
		  "steelblue": "#4682b4",
		  "tan": "#d2b48c",
		  "teal": "#008080",
		  "thistle": "#d8bfd8",
		  "tomato": "#ff6347",
		  "turquoise": "#40e0d0",
		  "violet": "#ee82ee",
		  "wheat": "#f5deb3",
		  "white": "#ffffff",
		  "whitesmoke": "#f5f5f5",
		  "aliceblue": "#f0f8ff",
		  "antiquewhite": "#faebd7",
		  "aqua": "#00ffff",
		  "aquamarine": "#7fffd4",
		  "azure": "#f0ffff",
		  "beige": "#f5f5dc",
		  "bisque": "#ffe4c4",
		  "blanchedalmond": "#ffebcd",
		  "blueviolet": "#8a2be2",
		  "yellowgreen": "#9acd32"
		};
		return colorObject;
	}
	
	//create drop down
	function createDropDown(FirstRow,tfdropdowndiv,dropdownid,TFs,onChange){
		var Keys=[FirstRow];
		TFs.sort();
		for (var tf of TFs){
			Keys.push(tf);
		}
		
		d3.select(tfdropdowndiv)
		.select("select")
		.remove();
		
		d3.select(tfdropdowndiv)
		.append("select")
		.attr("id",dropdownid)
		.on("change",onChange)
		.selectAll("option")
		.data(Keys)
		.enter()
		.append("option")
		.text(function(d){
			return d;
		})
		.attr("value",function(d){
			return d.__data__;
			});
	
	}
	
	function gedgeonchange(groupid){
		var edgeid=document.getElementById('g'+groupid+'edgeid').value
		var edgeid=edgeid.split(":")[0].slice(1);
		var gcolor1=document.getElementById("scolor1").value;
		var gcolor2=document.getElementById("scolor2").value;
		
		var se=d3.selectAll('.cells')[0]
		se=se.filter(function(d){
			var celledge=d.__data__.location[0];
			if (celledge==edgeid){
				return true;
			}return false;
		});
		
		//currentSelected=se;
		
		if (groupid==1){
			s1=s1.concat(se);
			s2=s2.filter(function(d){
				if (s1.indexOf(d)!=-1){
					return false;
				}return true;
			});
			console.log('test');

		}else{
			s2=s2.concat(se);
			s1=s1.filter(function(d){
				if (s2.indexOf(d)!=-1){
					return false;
				}return true;
			});
		}
		
		
		//currentSelected=currentSelected.concat(se);
		d3.selectAll(s1)
		.style("fill",gcolor1);
		document.getElementById("sgroup1").innerHTML=s1.length;
		
		d3.selectAll(s2)
		.style("fill",gcolor2);
		document.getElementById("sgroup2").innerHTML=s2.length;
		
		/*
		d3.selectAll(currentSelected)
		.style("fill",gcolor);
		document.getElementById("sgroup"+groupid).innerHTML=currentSelected.length;
		
		if (groupid==1){
			s1=currentSelected;
		}else{
			s2=currentSelected;
		}
		currentSelected=[];
		* */
		
	}
	
	
	// visualize all the cells in 2D space projected by t-sne and pca
	function vizCells(viztypep){
		//var plotType=document.getElementById(viztypep).checked;
		var plotType=$("input[name=viztype]:checked").val();
		var xData=[];
		var xLabels=[];
		var xTypeLabels=[];
		xcells=RL[2];
		
		for (var i in vizcells){
			var cell=vizcells[i];
			if (plotType=="pca"){
				xData.push(cell.PE);
			}else if (plotType=="dfmap"){
				xData.push(cell.ME)
			}else {
				xData.push(cell.TE);
			}
			
			var e=edges[xcells[i].location[0]];
			//xLabels.push(nodes[e.source_node].id+"->"+nodes[e.end_node].id);
			xLabels.push('P'+(nodes[e.end_node].id-1));
			xTypeLabels.push(cell.typeLabel);
		}
		
		var colorby=$("#ColorCode").val();
		if (colorby=='time'){
			scatterplot(xData,xTypeLabels,xTypeLabels);
		}else{
			scatterplot(xData,xLabels,xTypeLabels);
		}
	};
	
	function scatterplot(xData,xLabels,xTypeLabels){
			var dlabels={};
			for (var ix in xData){
				var di=xData[ix];
				var li=xLabels[ix];
				if (!(li in dlabels)){
					dlabels[li]=[{"x": di[0],"y": di[1]}];
				}else{
					dlabels[li].push({"x": di[0],"y": di[1]});
				}
			}
				
			var CombinedDataSet=[];
			var colorCounter=0;
			for (var di in dlabels){
				var xlabel=di;
				var xList=dlabels[di];
				CombinedDataSet.push({"label": xlabel, "data": xList, "backgroundColor": colorList[colorCounter]});
				colorCounter+=1;
			}
			
			var plotdata={"datasets":CombinedDataSet};
			
			var newW3 = open('','_blank','height=1000,width=1200')
			newW3.document.write('<head><title>scatter plot</title> </head><body></body>');
			d3.select(newW3.document.body).append("canvas")
						.attr("id","scatterplot")
						.attr("width","1000px")
						.attr("height","800px");
						
			var ctx=newW3.document.getElementById("scatterplot").getContext('2d');
			
			
			var scatterChart = new Chart(ctx, {
				type: 'scatter',
				data: plotdata,
				options: {
				responsive : false,
				scales: {
					xAxes: [{
						type: 'linear',
						position: 'bottom'
					}]
				},
				tooltips: {
					callbacks:{
						label: function(tooltipItem,data){
							var pos=[tooltipItem.xLabel,tooltipItem.yLabel];
							for (var ix in  xData){
								if (xData[ix].toString()==pos.toString()){
									return xTypeLabels[ix];
								}							
							}
						}
					}
				}
			}
				
			});	
		}
		
		//
		
		function panthergoInput(keys,species){
		var species=`
		<body>
		<span> choose your species: </span> <br>
		<select
		   id="rte_species"
		   class="form-control"
		   name="species">
		  <option value="HUMAN">Homo sapiens</option>
		  <option value="MOUSE">Mus musculus</option>
		  <option value="RAT">Rattus norvegicus</option>
		  <option value="CHICK">Gallus gallus</option>
		  <option value="DANRE">Danio rerio</option>
		  <option value="DROME">Drosophila melanogaster</option>
		  <option value="CAEEL">Caenorhabditis elegans</option>
		  <option value="YEAST">Saccharomyces cerevisiae</option>
		  <option value="SCHPO">Schizosaccharomyces pombe</option>
		  <option value="DICDI">Dictyostelium discoideum</option>
		  <option value="ARATH">Arabidopsis thaliana</option>
		  <option value="ECOLI">Escherichia coli</option>
		  <option value="EMENI">Emericella nidulans</option>
		  <option value="ANOCA">Anolis carolinensis</option>
		  <option value="ANOGA">Anopheles gambiae</option>
		  <option value="AQUAE">Aquifex aeolicus</option>
		  <option value="ASHGO">Ashbya gossypii</option>
		  <option value="BACCR">Bacillus cereus</option>
		  <option value="BACSU">Bacillus subtilis</option>
		  <option value="BACTN">Bacteroides thetaiotaomicron</option>
		  <option value="BATDJ">Batrachochytrium dendrobatidis</option>
		  <option value="BOVIN">Bos taurus</option>
		  <option value="BRADI">Brachypodium distachyon</option>
		  <option value="BRAJA">Bradyrhizobium japonicum</option>
		  <option value="BRAFL">Branchiostoma floridae</option>
		  <option value="CAEBR">Caenorhabditis briggsae</option>
		  <option value="CANAL">Candida albicans</option>
		  <option value="CANFA">Canis familiaris</option>
		  <option value="CHLTR">Chlamydia trachomatis</option>
		  <option value="CHLRE">Chlamydomonas reinhardtii</option>
		  <option value="CHLAA">Chloroflexus aurantiacus</option>
		  <option value="CIOIN">Ciona intestinalis</option>
		  <option value="CLOBH">Clostridium botulinum</option>
		  <option value="COXBU">Coxiella burnetii</option>
		  <option value="CRYNJ">Cryptococcus neoformans</option>
		  <option value="DAPPU">Daphnia pulex</option>
		  <option value="DEIRA">Deinococcus radiodurans</option>
		  <option value="DICTD">Dictyoglomus turgidum</option>
		  <option value="DICPU">Dictyostelium purpureum</option>
		  <option value="ENTHI">Entamoeba histolytica</option>
		  <option value="HORSE">Equus caballus</option>
		  <option value="FELCA">Felis catus</option>
		  <option value="FUSNN">Fusobacterium nucleatum</option>
		  <option value="GEOSL">Geobacter sulfurreducens</option>
		  <option value="GIAIC">Giardia intestinalis</option>
		  <option value="GLOVI">Gloeobacter violaceus</option>
		  <option value="SOYBN">Glycine max</option>
		  <option value="HAEIN">Haemophilus influenzae</option>
		  <option value="HALSA">Halobacterium salinarum</option>
		  <option value="IXOSC">Ixodes scapularis</option>
		  <option value="KORCO">Korarchaeum cryptofilum</option>
		  <option value="LEIMA">Leishmania major</option>
		  <option value="LEPIN">Leptospira interrogans</option>
		  <option value="LISMO">Listeria monocytogenes</option>
		  <option value="MACMU">Macaca mulatta</option>
		  <option value="METJA">Methanocaldococcus jannaschii</option>
		  <option value="METAC">Methanosarcina acetivorans</option>
		  <option value="MONDO">Monodelphis domestica</option>
		  <option value="MONBE">Monosiga brevicollis</option>
		  <option value="MYCTU">Mycobacterium tuberculosis</option>
		  <option value="NEMVE">Nematostella vectensis</option>
		  <option value="NEUCR">Neurospora crassa</option>
		  <option value="ORNAN">Ornithorhynchus anatinus</option>
		  <option value="ORYSJ">Oryza sativa</option>
		  <option value="PANTR">Pan troglodytes</option>
		  <option value="PHANO">Phaeosphaeria nodorum</option>
		  <option value="PHYPA">Physcomitrella patens</option>
		  <option value="PHYIT">Phytophthora infestans</option>
		  <option value="PLAF7">Plasmodium falciparum</option>
		  <option value="POPTR">Populus trichocarpa</option>
		  <option value="PRIPA">Pristionchus pacificus</option>
		  <option value="PSEAE">Pseudomonas aeruginosa</option>
		  <option value="PUCGT">Puccinia graminis</option>
		  <option value="PYRAE">Pyrobaculum aerophilum</option>
		  <option value="PYRKO">Pyrococcus kodakaraensis</option>
		  <option value="RHOBA">Rhodopirellula baltica</option>
		  <option value="SALTY">Salmonella typhimurium</option>
		  <option value="SCHMA">Schistosoma mansoni</option>
		  <option value="SCLS1">Sclerotinia sclerotiorum</option>
		  <option value="SHEON">Shewanella oneidensis</option>
		  <option value="SOLLC">Solanum lycopersicum</option>
		  <option value="SORBI">Sorghum bicolor</option>
		  <option value="STAA8">Staphylococcus aureus</option>
		  <option value="STRR6">Streptococcus pneumoniae</option>
		  <option value="STRCO">Streptomyces coelicolor</option>
		  <option value="STRPU">Strongylocentrotus purpuratus</option>
		  <option value="SULSO">Sulfolobus solfataricus</option>
		  <option value="PIG">Sus scrofa</option>
		  <option value="SYNY3">Synechocystis</option>
		  <option value="TAKRU">Takifugu rubripes</option>
		  <option value="TETTS">Tetrahymena thermophila</option>
		  <option value="THAPS">Thalassiosira pseudonana</option>
		  <option value="THEYD">Thermodesulfovibrio yellowstonii</option>
		  <option value="THEMA">Thermotoga maritima</option>
		  <option value="TRIVA">Trichomonas vaginalis</option>
		  <option value="TRIAD">Trichoplax adhaerens</option>
		  <option value="brucei">TRYB2 Trypanosoma brucei</option>
		  <option value="USTMA">Ustilago maydis</option>
		  <option value="VIBCH">Vibrio cholerae</option>
		  <option value="VITVI">Vitis vinifera</option>
		  <option value="XANCP">Xanthomonas campestris</option>
		  <option value="XENTR">Xenopus tropicalis</option>
		  <option value="YARLI">Yarrowia lipolytica</option>
		  <option value="YERPE">Yersinia pestis</option>
		</select>
		<input type="submit" id="gospecies">
		</body>
		`
		var ww=open("",'_blank','height=300,width=400,scrollbars=yes');
		ww.document.write(species);
		d3.select(ww.document.body).select("#gospecies")
		.on("click",function(){
			var ss=ww.document.getElementById("rte_species").value;
			panthergoInputSubmit(keys,ss);
		});		
	}
	
	function panthergoInputSubmit(keys,species){
		var link_pre="http://pantherdb.org/webservices/go/overrep.jsp?input=";
		if (species==undefined){
			var link_suffix="&species=MOUSE";
		}else{
			var link_suffix="&species="+species;
		}
		
		keys=keys.join("\n");
		key=encodeURIComponent(keys);
		var link=link_pre+keys+link_suffix;
		link=encodeURI(link);
		var ww=open(link,'_blank','height=600,width=800,left=1200,top=200,scrollbars=yes');
	}
