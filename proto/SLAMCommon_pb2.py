# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: SLAMCommon.proto
"""Generated protocol buffer code."""
from google.protobuf.internal import builder as _builder
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()




DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n\x10SLAMCommon.proto\x12\nslamcommon\"5\n\x07Version\x12\x18\n\x10protocol_version\x18\x01 \x01(\t\x12\x10\n\x08model_id\x18\x02 \x01(\t\".\n\x06Header\x12\x12\n\ntime_stamp\x18\x01 \x01(\x05\x12\x10\n\x08\x66rame_id\x18\x02 \x01(\x05\"5\n\x05Point\x12\t\n\x01x\x18\x01 \x01(\x01\x12\t\n\x01y\x18\x02 \x01(\x01\x12\t\n\x01z\x18\x03 \x01(\x01\x12\x0b\n\x03std\x18\x04 \x01(\x01\"\xbd\x02\n\x05\x42ox3D\x12#\n\x08lower_lt\x18\x01 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08lower_lb\x18\x02 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08lower_rb\x18\x03 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08lower_rt\x18\x04 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08upper_lt\x18\x05 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08upper_lb\x18\x06 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08upper_rb\x18\x07 \x01(\x0b\x32\x11.slamcommon.Point\x12#\n\x08upper_rt\x18\x08 \x01(\x0b\x32\x11.slamcommon.Point\x12\x0c\n\x04\x63onf\x18\t \x01(\x02\"u\n\x07Polygon\x12!\n\x06points\x18\x01 \x03(\x0b\x32\x11.slamcommon.Point\x12!\n\x06normal\x18\x02 \x01(\x0b\x32\x11.slamcommon.Point\x12$\n\tdirection\x18\x03 \x01(\x0b\x32\x11.slamcommon.Point\"g\n\nORBFeature\x12\x0e\n\x06octave\x18\x01 \x01(\x05\x12\t\n\x01x\x18\x02 \x01(\x02\x12\t\n\x01y\x18\x03 \x01(\x02\x12\x0c\n\x04size\x18\x04 \x01(\x02\x12\x10\n\x08respones\x18\x05 \x01(\x02\x12\x13\n\x0b\x64\x65scriptors\x18\x06 \x03(\x05\x42\"\n com.mx.changan.datacollect.protob\x06proto3')

_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, globals())
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'SLAMCommon_pb2', globals())
if _descriptor._USE_C_DESCRIPTORS == False:

  DESCRIPTOR._options = None
  DESCRIPTOR._serialized_options = b'\n com.mx.changan.datacollect.proto'
  _VERSION._serialized_start=32
  _VERSION._serialized_end=85
  _HEADER._serialized_start=87
  _HEADER._serialized_end=133
  _POINT._serialized_start=135
  _POINT._serialized_end=188
  _BOX3D._serialized_start=191
  _BOX3D._serialized_end=508
  _POLYGON._serialized_start=510
  _POLYGON._serialized_end=627
  _ORBFEATURE._serialized_start=629
  _ORBFEATURE._serialized_end=732
# @@protoc_insertion_point(module_scope)
